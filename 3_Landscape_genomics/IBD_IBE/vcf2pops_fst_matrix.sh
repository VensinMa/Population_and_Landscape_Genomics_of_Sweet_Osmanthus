#!/bin/bash

# ==============================================================================
# 配置区域 (请在此处修改文件路径)
# ==============================================================================
# 输入文件
VCF_FILE="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/gea_core_loci/Core_Adaptive_SNPs_15925.recode.vcf"
POP_INFO_FILE="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/gea_core_loci/202samples.pop"

# 输出文件名称
RESULT_CSV="Pops_fst_results.csv"
MATRIX_CSV="Pops_fst_matrix.csv"

# 并行核心数
THREADS=20

# ==============================================================================
# 第一阶段：准备工作与 FST 计算 (Bash + VCFtools)
# ==============================================================================

# 检查依赖
if ! command -v env_parallel &> /dev/null; then
    echo "Error: GNU parallel (env_parallel) 未安装或未在 PATH 中。"
    exit 1
fi

source $(which env_parallel.bash)

echo ">>> 开始初始化文件夹..."
mkdir -p pop_files
mkdir -p log_files

# 1. 根据总群体信息文件创建每个群体的 .pop 文件
echo ">>> 正在拆分群体文件..."
awk '{print $1 > "pop_files/"$2".pop"}' "$POP_INFO_FILE"

# 2. 生成两两组合
pop_files=(pop_files/*.pop)
combinations=()

echo ">>> 生成群体组合..."
for (( i=0; i<${#pop_files[@]}; i++ )); do
    for (( j=i+1; j<${#pop_files[@]}; j++ )); do
        combinations+=("${pop_files[i]}:${pop_files[j]}")
    done
done

echo ">>> 共生成 ${#combinations[@]} 对组合。"

# 3. 定义计算 FST 的函数
calculate_fst() {
    pair=$1
    IFS=":" read -r group1 group2 <<< "$pair"
    base_name1=$(basename "${group1%.*}")
    base_name2=$(basename "${group2%.*}")
    output_base="fst_${base_name1}_vs_${base_name2}"
    output_file="${output_base}"
    
    # 计算 FST (VCFtools)
    vcftools --vcf "$VCF_FILE" \
             --weir-fst-pop "$group1" \
             --weir-fst-pop "$group2" \
             --fst-window-size 100000 \
             --fst-window-step 10000 \
             --out "log_files/$output_file" >& "log_files/$output_file.log"
    
    OUTPUT=$(cat "log_files/$output_file.log")

    # 提取 FST 估计值
    MEAN_FST=$(echo "$OUTPUT" | grep "Weir and Cockerham mean Fst estimate" | awk '{print $NF}')
    WEIGHTED_FST=$(echo "$OUTPUT" | grep "Weir and Cockerham weighted Fst estimate" | awk '{print $NF}')

    # 处理可能的空值（如果计算失败）
    if [ -z "$MEAN_FST" ]; then MEAN_FST="NA"; fi
    if [ -z "$WEIGHTED_FST" ]; then WEIGHTED_FST="NA"; fi

    # 输出结果
    echo "${base_name1},${base_name2},${MEAN_FST},${WEIGHTED_FST}" >> "$RESULT_CSV_PATH"
    echo "Done: $base_name1 vs $base_name2 (Mean Fst: $MEAN_FST)"
}

export -f calculate_fst
export VCF_FILE
export RESULT_CSV_PATH="$PWD/$RESULT_CSV" # 传递绝对路径以防并行环境路径问题

# 4. 初始化结果文件
echo "Group1,Group2,Mean_FST,Weighted_FST" > "$RESULT_CSV"

# 5. 并行计算
echo ">>> 开始并行计算 FST (使用 $THREADS 个核心)..."
env_parallel -j "$THREADS" calculate_fst ::: "${combinations[@]}"

echo ">>> FST 计算完成。结果保存在: $RESULT_CSV"

# ==============================================================================
# 第二阶段：转换矩阵 (嵌入式 R 脚本)
# ==============================================================================

echo ">>> 开始使用 R 生成矩阵..."

# 将 Bash 变量传递给 R，通过 Here-Doc (<< 'EOF') 嵌入 R 代码
# $1 对应 input_file (RESULT_CSV), $2 对应 output_file (MATRIX_CSV)

Rscript - "$RESULT_CSV" "$MATRIX_CSV" << 'EOF'
args <- commandArgs(trailingOnly = TRUE)
input_file  <- args[1]
output_file <- args[2]

# 读取数据
data <- read.csv(input_file, header=TRUE)

# 简单的错误检查
if (nrow(data) == 0) {
    stop("错误: 输入的 CSV 文件为空，无法生成矩阵。")
}

# 提取唯一的组名
unique_groups <- unique(c(as.character(data$Group1), as.character(data$Group2)))
n_groups <- length(unique_groups)

cat(sprintf("检测到 %d 个群体，正在构建矩阵...\n", n_groups))

# 创建全零矩阵
fst_matrix <- matrix(0, nrow=n_groups, ncol=n_groups, 
                     dimnames=list(unique_groups, unique_groups))

# 填充矩阵
for (i in 1:nrow(data)) {
    row_name <- as.character(data[i, "Group1"])
    col_name <- as.character(data[i, "Group2"])
    
    # 获取 Mean_FST，如果是 NA 则设为 0 或其他标记
    val <- as.numeric(data[i, "Mean_FST"])
    if(is.na(val)) val <- 0 
    
    # 填充对称位置
    fst_matrix[row_name, col_name] <- val
    fst_matrix[col_name, row_name] <- val 
}

# 保存矩阵
write.csv(fst_matrix, output_file)
cat("成功: FST 矩阵已保存到文件：", output_file, "\n")
EOF

echo ">>> 流程全部结束！"
