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
STEP_SIZE=10000     # 10kb

# ================= 2. 环境准备 =================
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit

echo "正在准备样品列表..."
mkdir -p sample_lists/population
mkdir -p sample_lists/lineage

# --- A. 生成 [群体 Population] (第2列) 列表 ---
# 获取所有不重复的群体名
awk '{print $2}' "$POP_FILE" | sort | uniq > sample_lists/pop_names.txt

while read -r POP_NAME; do
    # 提取该群体的 SampleID (第1列)
    awk -v p="$POP_NAME" '$2 == p {print $1}' "$POP_FILE" > "sample_lists/population/${POP_NAME}.txt"
done < sample_lists/pop_names.txt

# --- B. 生成 [谱系 Lineage] (第3列) 列表 ---
# 获取所有不重复的谱系名 (假设文件是制表符分隔，处理空格问题)
# 如果文件是空格分隔且谱系名里有空格，需要特别小心。这里假设使用 \t 或标准空格分割
# 为了安全，我们对文件名中的空格进行替换
awk -F'\t' '{print $3}' "$POP_FILE" | sort | uniq > sample_lists/lineage_names.txt

while read -r LINEAGE_NAME; do
    # 1. 创建安全的文件名 (把空格换成下划线)
    SAFE_NAME=$(echo "$LINEAGE_NAME" | tr ' ' '_')
    
    # 2. 提取该谱系的 SampleID
    # 注意：这里匹配原始的 LINEAGE_NAME
    awk -F'\t' -v l="$LINEAGE_NAME" '$3 == l {print $1}' "$POP_FILE" > "sample_lists/lineage/${SAFE_NAME}.txt"
done < sample_lists/lineage_names.txt


# ================= 3. 定义计算函数 =================
# 参数: DatasetName(ALL/LD) VCFPath Level(Population/Lineage)
run_calc() {
    local D_NAME=$1
    local VCF_IN=$2
    local LEVEL=$3
    
    local LIST_DIR="sample_lists/${LEVEL,,}" # 转小写: population 或 lineage
    local OUT_DIR="${WORKDIR}/${D_NAME}/${LEVEL}"
    
    echo "----------------------------------------------------"
    echo "正在计算: [${D_NAME}] - [${LEVEL} Level]"
    echo "输入列表: ${LIST_DIR}"
    echo "输出目录: ${OUT_DIR}"
    
    mkdir -p "$OUT_DIR"

    # 遍历列表目录下的所有 .txt 文件
    for list_file in "${LIST_DIR}"/*.txt; do
        # 获取文件名作为前缀 (去掉了路径和.txt后缀)
        group_name=$(basename "$list_file" .txt)
        
        # 跳过空文件
        if [ ! -s "$list_file" ]; then continue; fi

        echo "  -> 处理分组: ${group_name}"

        # 1. 计算 Pi
        vcftools --gzvcf "$VCF_IN" \
                 --keep "$list_file" \
                 --window-pi "$WINDOW_SIZE" \
                # --window-pi-step "$STEP_SIZE" \
                 --out "${OUT_DIR}/${group_name}" 2> /dev/null

        # 2. 计算 Tajima's D
        vcftools --gzvcf "$VCF_IN" \
                 --keep "$list_file" \
                 --TajimaD "$WINDOW_SIZE" \
                 --out "${OUT_DIR}/${group_name}" 2> /dev/null
    done
}

# ================= 4. 开始批量运行 =================

# 任务 1: ALL_SNP 数据集
run_calc "ALL_SNP" "$VCF_ALL" "Population"
run_calc "ALL_SNP" "$VCF_ALL" "Lineage"

# 任务 2: LD_SNP 数据集
# (注意: LD prune 后的数据算 Pi 会被低估，算 TajimaD 也不准，但通常也会算一下作为参考)
run_calc "LD_SNP" "$VCF_LD" "Population"
run_calc "LD_SNP" "$VCF_LD" "Lineage"

echo "===================================================="
echo "所有计算完成！"
echo "结果保存在: ${WORKDIR}/[ALL_SNP|LD_SNP]/[Population|Lineage]/"
