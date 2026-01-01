#!/bin/bash

# ================= 1. 参数配置区域 =================
# 工作主目录
WORKDIR="/home/vensin/workspace/snpcalling_wild/12.population_genetics/Genetic_diversity"

# 输入文件路径
VCF_ALL="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz"
VCF_LD="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.recode.vcf.gz"
POP_FILE="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202samples.pop"

# 滑窗参数 (仅用于 Pi 和 Tajima's D)
WINDOW_SIZE=100000  # 100kb

# ================= 2. 核心数与输入检查 =================

# --- 2.1 检查输入文件 (Fail Fast) ---
echo "正在检查输入文件..."
MISSING=0
if [ ! -f "$VCF_ALL" ]; then echo "Error: VCF_ALL 文件不存在: $VCF_ALL"; MISSING=1; fi
if [ ! -f "$VCF_LD" ]; then echo "Error: VCF_LD 文件不存在: $VCF_LD"; MISSING=1; fi
if [ ! -f "$POP_FILE" ]; then echo "Error: POP_FILE 文件不存在: $POP_FILE"; MISSING=1; fi

if [ "$MISSING" -eq 1 ]; then
    echo "请检查文件路径配置，脚本已终止。"
    exit 1
fi
echo "输入文件检查通过。"

# --- 2.2 自动计算并行任务数 ---
TOTAL_CORES=$(nproc)
# 默认使用 80% 的核心数，防止卡死服务器
MAX_JOBS=$(( TOTAL_CORES * 8 / 10 ))
if [ "$MAX_JOBS" -lt 1 ]; then MAX_JOBS=1; fi

echo "检测到系统核心数: ${TOTAL_CORES}"
echo "并行任务数设置为: ${MAX_JOBS}"

# ================= 3. 环境准备 =================
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit

JOB_FILE="${WORKDIR}/pending_jobs.txt"
PROGRESS_LOG="${WORKDIR}/progress.log"
> "$JOB_FILE"
> "$PROGRESS_LOG"

echo "正在生成样品列表..."
mkdir -p sample_lists/population
mkdir -p sample_lists/lineage
mkdir -p sample_lists/species

# [优化] 统一使用 -F'\t' 确保精确解析
# --- A. Population (第2列) ---
awk -F'\t' '{print $2}' "$POP_FILE" | sort | uniq > sample_lists/pop_names.txt
while read -r POP_NAME; do
    # 移除特殊字符，确保文件名安全
    SAFE_POP=$(echo "$POP_NAME" | tr -d '[:space:]' | tr -cd '[:alnum:]_')
    awk -F'\t' -v p="$POP_NAME" '$2 == p {print $1}' "$POP_FILE" > "sample_lists/population/${POP_NAME}.txt"
done < sample_lists/pop_names.txt

# --- B. Lineage (第3列) ---
awk -F'\t' '{print $3}' "$POP_FILE" | sort | uniq > sample_lists/lineage_names.txt
while read -r LINEAGE_NAME; do
    SAFE_NAME=$(echo "$LINEAGE_NAME" | tr ' ' '_')
    awk -F'\t' -v l="$LINEAGE_NAME" '$3 == l {print $1}' "$POP_FILE" > "sample_lists/lineage/${SAFE_NAME}.txt"
done < sample_lists/lineage_names.txt

# --- C. Species (用于Pi/TD的全样本计算) ---
awk -F'\t' '{print $1}' "$POP_FILE" > "sample_lists/species/All_Samples.txt"

# ================= 4. 定义命令生成函数 =================

# --- 函数 A: 分组计算 (仅用于 Pi 和 Tajima's D) ---
# 逻辑: 这些指标依赖于群体划分，必须按 Pop/Lineage 循环运行
run_calc_gen() {
    local D_NAME=$1
    local VCF_IN=$2
    local LEVEL=$3
    
    local LIST_DIR="sample_lists/${LEVEL,,}" # 转小写
    local OUT_DIR="${WORKDIR}/${D_NAME}/${LEVEL}"
    
    mkdir -p "$OUT_DIR"

    for list_file in "${LIST_DIR}"/*.txt; do
        group_name=$(basename "$list_file" .txt)
        if [ ! -s "$list_file" ]; then continue; fi
        
        LOG_FILE="${OUT_DIR}/${group_name}.log"

        # 1. Pi (Nucleotide Diversity) - 滑窗
        CMD_PI="vcftools --gzvcf $VCF_IN --keep $list_file --window-pi $WINDOW_SIZE --out ${OUT_DIR}/${group_name} >> ${LOG_FILE} 2>&1"
        echo "${CMD_PI} && echo done >> ${PROGRESS_LOG}" >> "$JOB_FILE"
        
        # 2. Tajima's D - 滑窗
        CMD_TD="vcftools --gzvcf $VCF_IN --keep $list_file --TajimaD $WINDOW_SIZE --out ${OUT_DIR}/${group_name} >> ${LOG_FILE} 2>&1"
        echo "${CMD_TD} && echo done >> ${PROGRESS_LOG}" >> "$JOB_FILE"
        
        # 注意: 这里不再计算 het 和 hardy
    done
}

# --- 函数 B: 全局计算 (用于 Ho, He, Fis) ---
# 逻辑: 这些指标基于所有个体计算，不需要 --keep 参数，只运行一次
run_global_stats() {
    local D_NAME=$1
    local VCF_IN=$2
    
    # 结果存放在单独的 Global_Stats 目录
    local OUT_DIR="${WORKDIR}/${D_NAME}/Global_Stats"
    mkdir -p "$OUT_DIR"
    
    local LOG_FILE="${OUT_DIR}/global_calc.log"
    local OUT_PREFIX="${OUT_DIR}/All_Samples"

    # 3. Ho & He (Site-based) - 全局
    # 不加 --keep 默认使用 VCF 中所有样本
    CMD_HARDY="vcftools --gzvcf $VCF_IN --hardy --out ${OUT_PREFIX} >> ${LOG_FILE} 2>&1"
    echo "${CMD_HARDY} && echo done >> ${PROGRESS_LOG}" >> "$JOB_FILE"

    # 4. Ho & He (Individual-based) - 全局
    # 这样计算出的 Expected Homozygosity 是基于整个大群体的频率
    CMD_HET="vcftools --gzvcf $VCF_IN --het --out ${OUT_PREFIX} >> ${LOG_FILE} 2>&1"
    echo "${CMD_HET} && echo done >> ${PROGRESS_LOG}" >> "$JOB_FILE"
}

# ================= 5. 生成所有任务 =================
echo "正在生成任务列表..."

# --- 1. 分组任务 (Pi, Tajima's D) ---
run_calc_gen "ALL_SNP" "$VCF_ALL" "Population"
run_calc_gen "ALL_SNP" "$VCF_ALL" "Lineage"
run_calc_gen "ALL_SNP" "$VCF_ALL" "Species" # 这里 Species 组只算 Pi/TD

run_calc_gen "LD_SNP" "$VCF_LD" "Population"
run_calc_gen "LD_SNP" "$VCF_LD" "Lineage"
run_calc_gen "LD_SNP" "$VCF_LD" "Species"

# --- 2. 全局任务 (Het, Hardy) ---
# 针对两个数据集分别计算一次全局指标
echo "添加全局统计任务 (Het, Hardy)..."
run_global_stats "ALL_SNP" "$VCF_ALL"
run_global_stats "LD_SNP" "$VCF_LD"

TOTAL_TASKS=$(wc -l < "$JOB_FILE")
echo "任务生成完毕，共计 ${TOTAL_TASKS} 个任务。"

# ================= 6. 并行执行与监控 =================

monitor_progress() {
    local total=$1
    local start_time=$(date +%s)
    # 尝试隐藏光标 (如果终端支持)
    tput civis 2>/dev/null
    
    while true; do
        if [ -f "$PROGRESS_LOG" ]; then completed=$(wc -l < "$PROGRESS_LOG"); else completed=0; fi
        if [ "$total" -gt 0 ]; then percent=$(( completed * 100 / total )); else percent=0; fi
        current_time=$(date +%s)
        elapsed=$(( current_time - start_time ))
        printf "\rProgress: [ %d / %d ] %d%% (Time: %ds) " "$completed" "$total" "$percent" "$elapsed"
        if [ "$completed" -ge "$total" ]; then break; fi
        sleep 1
    done
    
    # 恢复光标
    tput cnorm 2>/dev/null
    echo ""
}

echo "===================================================="
echo "开始并行处理 (并发数: ${MAX_JOBS})"
echo "===================================================="

monitor_progress "$TOTAL_TASKS" &
MONITOR_PID=$!

cat "$JOB_FILE" | xargs -P "$MAX_JOBS" -I {} sh -c "{}"

wait $MONITOR_PID

# ================= 7. 清理与完成 =================
rm "$JOB_FILE"
rm "$PROGRESS_LOG"

echo "===================================================="
echo "所有计算完成！"
echo "输出结构说明:"
echo "1. Pi/TajimaD 位于: [Dataset]/[Level]/ (按分组)"
echo "2. Ho/He/Fis  位于: [Dataset]/Global_Stats/ (基于全样本)"
echo "   - .het 文件包含了每个个体基于全群体频率的 F 值"
