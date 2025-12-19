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
STEP_SIZE=10000     # 10kb (如需启用步长，请修改下方命令生成部分)

# ================= 2. 核心数与输入检查 =================

# --- 2.1 检查输入文件是否存在 (Fail Fast) ---
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
# 获取系统逻辑核心数
TOTAL_CORES=$(nproc)
# 默认使用 80% 的核心数
MAX_JOBS=$(( TOTAL_CORES * 8 / 10 ))
# 如果计算结果小于1，默认为1
if [ "$MAX_JOBS" -lt 1 ]; then MAX_JOBS=1; fi

# 【可选】如果你想手动指定，取消下面这行的注释并修改数字：
# MAX_JOBS=10

echo "检测到系统核心数: ${TOTAL_CORES}"
echo "并行任务数设置为: ${MAX_JOBS}"

# ================= 3. 环境准备 =================
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit

# 定义文件路径
JOB_FILE="${WORKDIR}/pending_jobs.txt"    # 存放待执行命令
PROGRESS_LOG="${WORKDIR}/progress.log"    # 存放进度标记
> "$JOB_FILE"      # 清空任务文件
> "$PROGRESS_LOG"  # 清空进度文件

echo "正在生成样品列表..."
mkdir -p sample_lists/population
mkdir -p sample_lists/lineage

# --- A. 生成 Population 列表 ---
awk '{print $2}' "$POP_FILE" | sort | uniq > sample_lists/pop_names.txt
while read -r POP_NAME; do
    awk -v p="$POP_NAME" '$2 == p {print $1}' "$POP_FILE" > "sample_lists/population/${POP_NAME}.txt"
done < sample_lists/pop_names.txt

# --- B. 生成 Lineage 列表 ---
awk -F'\t' '{print $3}' "$POP_FILE" | sort | uniq > sample_lists/lineage_names.txt
while read -r LINEAGE_NAME; do
    SAFE_NAME=$(echo "$LINEAGE_NAME" | tr ' ' '_')
    awk -F'\t' -v l="$LINEAGE_NAME" '$3 == l {print $1}' "$POP_FILE" > "sample_lists/lineage/${SAFE_NAME}.txt"
done < sample_lists/lineage_names.txt

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

        # 生成命令逻辑：
        # 1. 执行计算
        # 2. 成功后 (&&) 向 PROGRESS_LOG 文件追加一行 "done"
        # 3. 输出重定向到 /dev/null 防止屏幕刷屏
        
        # Pi 命令
        CMD_PI="vcftools --gzvcf $VCF_IN --keep $list_file --window-pi $WINDOW_SIZE --out ${OUT_DIR}/${group_name} > /dev/null 2>&1"
        echo "${CMD_PI} && echo done >> ${PROGRESS_LOG}" >> "$JOB_FILE"
        
        # Tajima's D 命令
        CMD_TD="vcftools --gzvcf $VCF_IN --keep $list_file --TajimaD $WINDOW_SIZE --out ${OUT_DIR}/${group_name} > /dev/null 2>&1"
        echo "${CMD_TD} && echo done >> ${PROGRESS_LOG}" >> "$JOB_FILE"
    done
}

# ================= 5. 生成所有任务 =================
echo "正在生成任务列表..."
run_calc_gen "ALL_SNP" "$VCF_ALL" "Population"
run_calc_gen "ALL_SNP" "$VCF_ALL" "Lineage"
run_calc_gen "LD_SNP" "$VCF_LD" "Population"
run_calc_gen "LD_SNP" "$VCF_LD" "Lineage"

TOTAL_TASKS=$(wc -l < "$JOB_FILE")
echo "任务生成完毕，共计 ${TOTAL_TASKS} 个任务。"

# ================= 6. 并行执行与进度监控 =================

# --- 定义进度条函数 ---
monitor_progress() {
    local total=$1
    local start_time=$(date +%s)
    
    while true; do
        # 计算已完成行数
        if [ -f "$PROGRESS_LOG" ]; then
            completed=$(wc -l < "$PROGRESS_LOG")
        else
            completed=0
        fi

        # 计算百分比
        if [ "$total" -gt 0 ]; then
            percent=$(( completed * 100 / total ))
        else
            percent=0
        fi
        
        # 计算耗时
        current_time=$(date +%s)
        elapsed=$(( current_time - start_time ))
        
        # 打印进度条 (使用 \r 回车不换行，实现原地刷新)
        # 格式: Progress: [Completed/Total] Percentage% (Elapsed Time)
        printf "\rProgress: [ %d / %d ] %d%% (Time: %ds) " "$completed" "$total" "$percent" "$elapsed"
        
        # 如果完成数等于总数，退出循环
        if [ "$completed" -ge "$total" ]; then
            break
        fi
        
        # 每 0.5 秒刷新一次
        sleep 0.5
    done
    echo "" # 换行
}

echo "===================================================="
echo "开始并行处理 (并发数: ${MAX_JOBS})"
echo "===================================================="

# 1. 在后台启动进度监控
monitor_progress "$TOTAL_TASKS" &
MONITOR_PID=$!

# 2. 开始执行并行任务
# xargs 说明: -P 并行数, -I 占位符, sh -c 执行具体的命令字符串
cat "$JOB_FILE" | xargs -P "$MAX_JOBS" -I {} sh -c "{}"

# 3. 等待监控进程结束 (确保进度条走到100%)
wait $MONITOR_PID

# ================= 7. 清理与完成 =================
rm "$JOB_FILE"
rm "$PROGRESS_LOG"

echo "===================================================="
echo "所有计算完成！"
echo "结果保存在: ${WORKDIR}"
