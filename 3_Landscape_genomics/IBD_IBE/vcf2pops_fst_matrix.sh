#!/bin/bash

# ==============================================================================
# 配置区域
# ==============================================================================

# 1. 输入文件路径
# 注意：脚本会自动检测文件后缀。如果是 .gz，会自动使用 --gzvcf
VCF_NEUTRAL="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/208_samples_snp_filtered.LD.pruned.recode.vcf.gz"
VCF_ADAPTIVE="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/gea_core_loci/Core_Adaptive_SNPs_15925.recode.vcf"

# 2. 群体信息文件
POP_INFO_FILE="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/gea_core_loci/202samples.pop"

# 3. 并行设置
THREADS=20

# ==============================================================================
# 核心函数
# ==============================================================================

# 检查依赖
if ! command -v env_parallel &> /dev/null; then
    echo "Error: GNU parallel (env_parallel) 未安装。"
    exit 1
fi
source $(which env_parallel.bash)

# FST 计算单元 (输出到 Stdout，不直接写文件)
calculate_fst_unit() {
    pair=$1
    # 环境变量: $CURRENT_VCF, $CURRENT_LOG_DIR, $VCF_FLAG (区分 --vcf/--gzvcf)
    
    IFS=":" read -r group1 group2 <<< "$pair"
    base_name1=$(basename "${group1%.*}")
    base_name2=$(basename "${group2%.*}")
    
    # 定义临时日志路径
    output_base="fst_${base_name1}_vs_${base_name2}"
    log_file="${CURRENT_LOG_DIR}/${output_base}"
    
    # 运行 VCFtools
    # 错误日志重定向到 /dev/null 以保持屏幕整洁，除非你想调试
    vcftools "$VCF_FLAG" "$CURRENT_VCF" \
             --weir-fst-pop "$group1" \
             --weir-fst-pop "$group2" \
             --fst-window-size 100000 \
             --fst-window-step 10000 \
             --out "$log_file" > "${log_file}.console.log" 2>&1
    
    # 检查结果
    if [ -f "${log_file}.log" ]; then
        OUTPUT=$(cat "${log_file}.log")
        MEAN_FST=$(echo "$OUTPUT" | grep "Weir and Cockerham mean Fst estimate" | awk '{print $NF}')
        WEIGHTED_FST=$(echo "$OUTPUT" | grep "Weir and Cockerham weighted Fst estimate" | awk '{print $NF}')
        
        # 成功提取后，清理临时文件 (节省空间)
        rm "${log_file}.log" "${log_file}.console.log"
    else
        MEAN_FST="NA"
        WEIGHTED_FST="NA"
        # 失败时保留 console log 以便排查
        echo "Warning: Calculation failed for $base_name1 vs $base_name2" >&2
    fi

    # 处理空值
    if [ -z "$MEAN_FST" ]; then MEAN_FST="NA"; fi
    if [ -z "$WEIGHTED_FST" ]; then WEIGHTED_FST="NA"; fi

    # 输出CSV行 (标准输出)
    echo "${base_name1},${base_name2},${MEAN_FST},${WEIGHTED_FST}"
}
export -f calculate_fst_unit

# R 转换函数
run_r_matrix_conversion() {
    local in_csv=$1
    local out_csv=$2
    echo ">>> [R] 转换矩阵: $out_csv"
    Rscript - "$in_csv" "$out_csv" << 'EOF'
args <- commandArgs(trailingOnly = TRUE)
data <- read.csv(args[1], header=TRUE)
output_file <- args[2]

if (nrow(data) == 0) stop("CSV is empty.")

groups <- unique(c(as.character(data$Group1), as.character(data$Group2)))
mat <- matrix(0, nrow=length(groups), ncol=length(groups), dimnames=list(groups, groups))

for (i in 1:nrow(data)) {
    r <- as.character(data[i, "Group1"]); c <- as.character(data[i, "Group2"])
    val <- as.numeric(data[i, "Mean_FST"])
    if(is.na(val)) val <- 0 
    mat[r, c] <- val; mat[c, r] <- val 
}
write.csv(mat, output_file)
EOF
}

# ==============================================================================
# 主流程
# ==============================================================================

# 1. 准备群体文件
echo ">>> [Init] 初始化群体文件..."
mkdir -p pop_files
awk '{print $1 > "pop_files/"$2".pop"}' "$POP_INFO_FILE"

# 生成组合
pop_files=(pop_files/*.pop)
combinations=()
for (( i=0; i<${#pop_files[@]}; i++ )); do
    for (( j=i+1; j<${#pop_files[@]}; j++ )); do
        combinations+=("${pop_files[i]}:${pop_files[j]}")
    done
done
echo ">>> [Init] 生成 ${#combinations[@]} 对组合。"

# 2. 管道函数
run_pipeline() {
    local TYPE=$1
    local VCF=$2
    
    echo -e "\n========================================================"
    echo ">>> 任务开始: $TYPE"
    echo ">>> 文件: $VCF"
    echo "========================================================"

    if [ ! -f "$VCF" ]; then echo "Error: 文件不存在 $VCF"; return; fi

    local LOG_DIR="log_files/${TYPE}"
    local RES_CSV="Pops_fst_results_${TYPE}.csv"
    local MAT_CSV="Pops_fst_matrix_${TYPE}.csv"

    mkdir -p "$LOG_DIR"
    
    # 写入 CSV 标题
    echo "Group1,Group2,Mean_FST,Weighted_FST" > "$RES_CSV"

    # 设置 VCF 标志 (自动检测 .gz)
    export CURRENT_VCF="$VCF"
    export CURRENT_LOG_DIR="$LOG_DIR"
    if [[ "$VCF" == *.gz ]]; then
        export VCF_FLAG="--gzvcf"
        echo ">>> 检测到压缩文件，使用 --gzvcf 模式"
    else
        export VCF_FLAG="--vcf"
    fi

    # 并行计算 (输出直接追加到 CSV，注意 --bar 显示进度条)
    echo ">>> [Parallel] 计算中..."
    env_parallel --bar -j "$THREADS" calculate_fst_unit ::: "${combinations[@]}" >> "$RES_CSV"

    # R 转换
    run_r_matrix_conversion "$RES_CSV" "$MAT_CSV"
    echo ">>> $TYPE 完成。结果: $MAT_CSV"
}

# ==============================================================================
# 执行
# ==============================================================================

run_pipeline "Neutral" "$VCF_NEUTRAL"
run_pipeline "Adaptive" "$VCF_ADAPTIVE"

echo -e "\n>>> 全部完成！"
