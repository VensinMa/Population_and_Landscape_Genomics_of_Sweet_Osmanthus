#!/bin/bash

# ================= 1. 参数配置区域 =================
# 工作主目录
WORKDIR="/home/vensin/workspace/snpcalling_wild/12.population_genetics/Genetic_diversity"

# 输入文件路径
VCF_ALL="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz"
VCF_LD="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.recode.vcf.gz"
POP_FILE="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202samples.pop"

# 滑窗参数
WINDOW_SIZE=100000  # 100kb
STEP_SIZE=10000     # 10kb (如果需要启用步长，请在下方生成命令处取消注释)

# 并行设置
# 获取系统逻辑核心数
TOTAL_CORES=$(nproc)
# 计算 80% 的核心数 (向下取整)
MAX_JOBS=$(( TOTAL_CORES * 8 / 10 ))
# 如果计算结果小于1，则默认为1
if [ "$MAX_JOBS" -lt 1 ]; then MAX_JOBS=1; fi

# 也可以在这里直接指定并行数，例如：
MAX_JOBS=24

echo "检测到系统核心数: ${TOTAL_CORES}"
echo "并行任务数设置为: ${MAX_JOBS} (占用约80%资源)"

# ================= 2. 输入文件检查 & 环境准备 =================
# 检查输入文件是否存在
if [ ! -f "$VCF_ALL" ]; then echo "Error: VCF_ALL 文件不存在: $VCF_ALL"; exit 1; fi
if [ ! -f "$VCF_LD" ]; then echo "Error: VCF_LD 文件不存在: $VCF_LD"; exit 1; fi
if [ ! -f "$POP_FILE" ]; then echo "Error: POP_FILE 文件不存在: $POP_FILE"; exit 1; fi

mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit

# 定义任务列表文件 (用于存放待执行的命令)
JOB_FILE="${WORKDIR}/pending_jobs.txt"
> "$JOB_FILE" # 清空或创建任务文件

echo "正在准备样品列表..."
mkdir -p sample_lists/population
mkdir -p sample_lists/lineage

# --- A. 生成 [群体 Population] (第2列) 列表 ---
awk '{print $2}' "$POP_FILE" | sort | uniq > sample_lists/pop_names.txt
while read -r POP_NAME; do
    awk -v p="$POP_NAME" '$2 == p {print $1}' "$POP_FILE" > "sample_lists/population/${POP_NAME}.txt"
done < sample_lists/pop_names.txt

# --- B. 生成 [谱系 Lineage] (第3列) 列表 ---
awk -F'\t' '{print $3}' "$POP_FILE" | sort | uniq > sample_lists/lineage_names.txt
while read -r LINEAGE_NAME; do
    SAFE_NAME=$(echo "$LINEAGE_NAME" | tr ' ' '_')
    # 注意：这里匹配原始的 LINEAGE_NAME，输出到 SAFE_NAME
    awk -F'\t' -v l="$LINEAGE_NAME" '$3 == l {print $1}' "$POP_FILE" > "sample_lists/lineage/${SAFE_NAME}.txt"
done < sample_lists/lineage_names.txt


# ================= 3. 定义命令生成函数 =================
# 注意：现在这个函数不再直接运行 vcftools，而是生成命令到 JOB_FILE
run_calc_gen() {
    local D_NAME=$1
    local VCF_IN=$2
    local LEVEL=$3
    
    local LIST_DIR="sample_lists/${LEVEL,,}" # 转小写
    local OUT_DIR="${WORKDIR}/${D_NAME}/${LEVEL}"
    
    echo "----------------------------------------------------"
    echo "正在生成任务: [${D_NAME}] - [${LEVEL} Level]"
    
    mkdir -p "$OUT_DIR"

    # 遍历列表目录下的所有 .txt 文件
    for list_file in "${LIST_DIR}"/*.txt; do
        group_name=$(basename "$list_file" .txt)
        if [ ! -s "$list_file" ]; then continue; fi

        # 生成 Pi 命令
        # 注意：这里把命令 echo 到文件里
        # 1. Pi 计算命令
        echo "vcftools --gzvcf $VCF_IN --keep $list_file --window-pi $WINDOW_SIZE --out ${OUT_DIR}/${group_name} > /dev/null 2>&1" >> "$JOB_FILE"
        
        # 2. Tajima's D 计算命令
        echo "vcftools --gzvcf $VCF_IN --keep $list_file --TajimaD $WINDOW_SIZE --out ${OUT_DIR}/${group_name} > /dev/null 2>&1" >> "$JOB_FILE"
    done
}

# ================= 4. 生成所有任务 =================

# 生成 ALL_SNP 数据集任务
run_calc_gen "ALL_SNP" "$VCF_ALL" "Population"
run_calc_gen "ALL_SNP" "$VCF_ALL" "Lineage"

# 生成 LD_SNP 数据集任务
run_calc_gen "LD_SNP" "$VCF_LD" "Population"
run_calc_gen "LD_SNP" "$VCF_LD" "Lineage"

# ================= 5. 并行执行任务 =================

TOTAL_JOBS=$(wc -l < "$JOB_FILE")
echo "===================================================="
echo "总共生成了 ${TOTAL_JOBS} 个 vcftools 任务。"
echo "开始并行执行 (并发数: ${MAX_JOBS})..."
echo "请稍候..."

# 使用 xargs 进行并行处理
# -P: 指定最大进程数
# -I {}: 占位符
# sh -c "{}": 执行命令字符串
cat "$JOB_FILE" | xargs -P "$MAX_JOBS" -I {} sh -c "{}"

# 清理任务文件
rm "$JOB_FILE"

echo "===================================================="
echo "所有计算完成！"
echo "结果保存在: ${WORKDIR}/[ALL_SNP|LD_SNP]/[Population|Lineage]/"
