#!/bin/bash

# ======================================================================
# GATK SelectVariants 并行分离 SNP 和 INDEL 脚本 (分类存储版)
# 资源需求: 约 100G 内存 (50G x 2)
# ======================================================================

# --- 1. 目录和文件配置 ---
input_dir="/home/data/8.gatk_mergevcfs"
input_vcf="$input_dir/merged.all_chromosomes.vcf.gz"
output_dir="/home/data/9.gatk_selectvariants"

# 定义 SNP 和 INDEL 的独立子目录
snp_dir="$output_dir/SNP"
indel_dir="$output_dir/INDEL"

reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
gatk_path="/home/vensin/software/gatk-4.6.2.0/gatk"

# 内存设置：每个任务 50G
JAVA_OPTS="-Xmx50g -Xms10g"

# --- 2. 准备工作 ---
# 创建主目录和子目录
mkdir -p "$output_dir" "$snp_dir" "$indel_dir"

# 主流程日志保存在最外层，方便查看总体进度
main_log="$output_dir/main_process.log"

# 定义子任务日志路径 (分别存放在各自的子文件夹中)
snp_log="$snp_dir/select_snp.log"
indel_log="$indel_dir/select_indel.log"

echo "=== GATK SelectVariants 并行任务开始: $(date) ===" | tee "$main_log"
echo "内存配置: 每个任务分配 $JAVA_OPTS (总计需 ~100G)" | tee -a "$main_log"
echo "SNP  输出目录: $snp_dir" | tee -a "$main_log"
echo "INDEL 输出目录: $indel_dir" | tee -a "$main_log"

# 检查输入文件
if [ ! -f "$input_vcf" ]; then
    echo "错误: 输入文件不存在: $input_vcf" | tee -a "$main_log"
    exit 1
fi

# --- 3. 启动 SNP 提取任务 (后台运行) ---
output_snp="$snp_dir/raw_snps.vcf.gz"

echo ">>> [后台任务 1] 启动 SNP 提取..." | tee -a "$main_log"
echo "    日志: $snp_log" | tee -a "$main_log"

(
    echo "开始提取 SNP: $(date)" > "$snp_log"
    $gatk_path --java-options "$JAVA_OPTS" SelectVariants \
        -R "$reference_genome" \
        -V "$input_vcf" \
        -select-type SNP \
        -O "$output_snp" \
        >> "$snp_log" 2>&1
    
    snp_exit_code=$?
    if [ $snp_exit_code -eq 0 ]; then
        echo "SNP 提取完成: $(date)" >> "$snp_log"
        exit 0
    else
        echo "SNP 提取失败!" >> "$snp_log"
        exit 1
    fi
) &
pid_snp=$!  # 记录 SNP 任务的进程 ID

# --- 4. 启动 INDEL 提取任务 (后台运行) ---
output_indel="$indel_dir/raw_indels.vcf.gz"

echo ">>> [后台任务 2] 启动 INDEL 提取..." | tee -a "$main_log"
echo "    日志: $indel_log" | tee -a "$main_log"

(
    echo "开始提取 INDEL: $(date)" > "$indel_log"
    $gatk_path --java-options "$JAVA_OPTS" SelectVariants \
        -R "$reference_genome" \
        -V "$input_vcf" \
        -select-type INDEL \
        -O "$output_indel" \
        >> "$indel_log" 2>&1
    
    indel_exit_code=$?
    if [ $indel_exit_code -eq 0 ]; then
        echo "INDEL 提取完成: $(date)" >> "$indel_log"
        exit 0
    else
        echo "INDEL 提取失败!" >> "$indel_log"
        exit 1
    fi
) &
pid_indel=$!  # 记录 INDEL 任务的进程 ID

# --- 5. 等待所有任务结束 ---
echo "" | tee -a "$main_log"
echo ">>> 所有任务已提交，正在后台运行..." | tee -a "$main_log"
echo "    SNP 任务 PID: $pid_snp" | tee -a "$main_log"
echo "    INDEL 任务 PID: $pid_indel" | tee -a "$main_log"
echo ">>> 请勿关闭终端，等待任务完成..." | tee -a "$main_log"

# 等待 SNP 任务
wait $pid_snp
status_snp=$?

# 等待 INDEL 任务
wait $pid_indel
status_indel=$?

echo "" | tee -a "$main_log"
echo "=== 任务执行报告 ===" | tee -a "$main_log"

# --- 6. 检查结果 ---
final_status=0

# 检查 SNP 结果
if [ $status_snp -eq 0 ] && [ -f "$output_snp" ]; then
    snp_size=$(ls -lh "$output_snp" | awk '{print $5}')
    echo "✓ SNP 提取成功" | tee -a "$main_log"
    echo "  文件: $output_snp" | tee -a "$main_log"
    echo "  大小: $snp_size" | tee -a "$main_log"
else
    echo "✗ SNP 提取失败 (退出码: $status_snp)" | tee -a "$main_log"
    echo "  请检查日志: $snp_log" | tee -a "$main_log"
    final_status=1
fi

# 检查 INDEL 结果
if [ $status_indel -eq 0 ] && [ -f "$output_indel" ]; then
    indel_size=$(ls -lh "$output_indel" | awk '{print $5}')
    echo "✓ INDEL 提取成功" | tee -a "$main_log"
    echo "  文件: $output_indel" | tee -a "$main_log"
    echo "  大小: $indel_size" | tee -a "$main_log"
else
    echo "✗ INDEL 提取失败 (退出码: $status_indel)" | tee -a "$main_log"
    echo "  请检查日志: $indel_log" | tee -a "$main_log"
    final_status=1
fi

# --- 7. 结束 ---
if [ $final_status -eq 0 ]; then
    echo "" | tee -a "$main_log"
    echo "全部完成! 下一步建议:" | tee -a "$main_log"
    echo "1. 进入 $snp_dir 对 SNP 进行硬过滤" | tee -a "$main_log"
    echo "2. 进入 $indel_dir 对 INDEL 进行硬过滤" | tee -a "$main_log"
    exit 0
else
    echo "" | tee -a "$main_log"
    echo "部分任务失败，请检查详细日志。" | tee -a "$main_log"
    exit 1
fi
