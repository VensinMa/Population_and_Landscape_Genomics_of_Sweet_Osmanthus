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
# 默认使用 80% 的核心数
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

# --- A. Population ---
awk '{print $2}' "$POP_FILE" | sort | uniq > sample_lists/pop_names.txt
while read -r POP_NAME; do
    awk -v p="$POP_NAME" '$2 == p {print $1}' "$POP_FILE" > "sample_lists/population/${POP_NAME}.txt"
done < sample_lists/pop_names.txt

# --- B. Lineage ---
awk -F'\t' '{print $3}' "$POP_FILE" | sort | uniq > sample_lists/lineage_names.txt
while read -r LINEAGE_NAME; do
    SAFE_NAME=$(echo "$LINEAGE_NAME" | tr ' ' '_')
    awk -F'\t' -v l="$LINEAGE_NAME" '$3 == l {print $1}' "$POP_FILE" > "sample_lists/lineage/${SAFE_NAME}.txt"
done < sample_lists/lineage_names.txt

# --- C. Species ---
awk '{print $1}' "$POP_FILE" > "sample_lists/species/All_Samples.txt"

# ================= 4. 定义命令生成函数 =================
run_calc_gen() {
    local D_NAME=$1
    local VCF_IN=$2
    local LEVEL=$3
    
    local LIST_DIR="sample_lists/${LEVEL,,}"
    local OUT_DIR="${WORKDIR}/${D_NAME}/${LEVEL}"
    
    mkdir -p "$OUT_DIR"

    for list_file in "${LIST_DIR}"/*.txt; do
        group_name=$(basename "$list_file" .txt)
        if [ ! -s "$list_file" ]; then continue; fi

        # 1. Pi (Nucleotide Diversity) - 滑窗
        CMD_PI="vcftools --gzvcf $VCF_IN --keep $list_file --window-pi $WINDOW_SIZE --out ${OUT_DIR}/${group_name} > /dev/null 2>&1"
        echo "${CMD_PI} && echo done >> ${PROGRESS_LOG}" >> "$JOB_FILE"
        
        # 2. Tajima's D - 滑窗
        CMD_TD="vcftools --gzvcf $VCF_IN --keep $list_file --TajimaD $WINDOW_SIZE --out ${OUT_DIR}/${group_name} > /dev/null 2>&1"
        echo "${CMD_TD} && echo done >> ${PROGRESS_LOG}" >> "$JOB_FILE"

        # 3. Ho & He (Site-based) - 按位点
        # 参数: --hardy
        # 输出后缀: .hwe
        # 用途: 通过每个位点的基因型计数计算 Ho/He
        CMD_HARDY="vcftools --gzvcf $VCF_IN --keep $list_file --hardy --out ${OUT_DIR}/${group_name} > /dev/null 2>&1"
        echo "${CMD_HARDY} && echo done >> ${PROGRESS_LOG}" >> "$JOB_FILE"

        # 4. Ho & He (Individual-based) - 按个体
        # 参数: --het
        # 输出后缀: .het
        # 用途: 通过每个个体的纯合子计数(O_HOM, E_HOM)计算 Ho/He
        CMD_HET="vcftools --gzvcf $VCF_IN --keep $list_file --het --out ${OUT_DIR}/${group_name} > /dev/null 2>&1"
        echo "${CMD_HET} && echo done >> ${PROGRESS_LOG}" >> "$JOB_FILE"

    done
}

# ================= 5. 生成所有任务 =================
echo "正在生成任务列表..."

run_calc_gen "ALL_SNP" "$VCF_ALL" "Population"
run_calc_gen "ALL_SNP" "$VCF_ALL" "Lineage"
run_calc_gen "ALL_SNP" "$VCF_ALL" "Species"

run_calc_gen "LD_SNP" "$VCF_LD" "Population"
run_calc_gen "LD_SNP" "$VCF_LD" "Lineage"
run_calc_gen "LD_SNP" "$VCF_LD" "Species"

TOTAL_TASKS=$(wc -l < "$JOB_FILE")
echo "任务生成完毕，共计 ${TOTAL_TASKS} 个任务。"

# ================= 6. 并行执行与监控 =================

monitor_progress() {
    local total=$1
    local start_time=$(date +%s)
    while true; do
        if [ -f "$PROGRESS_LOG" ]; then completed=$(wc -l < "$PROGRESS_LOG"); else completed=0; fi
        if [ "$total" -gt 0 ]; then percent=$(( completed * 100 / total )); else percent=0; fi
        current_time=$(date +%s)
        elapsed=$(( current_time - start_time ))
        printf "\rProgress: [ %d / %d ] %d%% (Time: %ds) " "$completed" "$total" "$percent" "$elapsed"
        if [ "$completed" -ge "$total" ]; then break; fi
        sleep 0.5
    done
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
echo "结果目录示例: ${WORKDIR}/ALL_SNP/Population/"
echo "文件说明:"
echo "  .windowed.pi -> Pi 值"
echo "  .Tajima.D -> Tajima's D 值"
echo "  .hwe -> Site-based Ho/He (源自 --hardy)"
echo "  .het -> Individual-based Ho/He (源自 --het)"
