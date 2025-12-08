#!/bin/bash

# 设置目录和文件
input_gvcf_dir="/home/vensin/workspace/snpcalling_wild/5.gatk_haplotypecaller/gvcf"
output_dir="/home/data/6.gatk_combinegvcfs"
reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
gatk_path="/home/vensin/software/gatk-4.6.2.0/gatk"

# 创建输出目录
mkdir -p "$output_dir"

# 设置日志文件
log_file="$output_dir/combinegvcfs.log"
echo "GATK CombineGVCFs 开始时间: $(date)" > "$log_file"

# 检查输入GVCF文件
gvcf_files=($(find "$input_gvcf_dir" -name '*.g.vcf.gz'))
gvcf_count=${#gvcf_files[@]}

if [ "$gvcf_count" -eq 0 ]; then
    echo "错误: 未找到GVCF文件" | tee -a "$log_file"
    exit 1
fi

echo "找到 $gvcf_count 个GVCF文件" >> "$log_file"

# 创建GVCF文件列表
gvcf_list_file="$output_dir/gvcf_files.list"
printf "%s\n" "${gvcf_files[@]}" | sort > "$gvcf_list_file"

# 设置合并参数
java_memory="100g"
combined_gvcf="$output_dir/combined.g.vcf.gz"

echo "开始合并 $gvcf_count 个GVCF文件..." >> "$log_file"

# 执行CombineGVCFs
$gatk_path --java-options "-Xmx${java_memory} -Xms80g" CombineGVCFs \
    -R "$reference_genome" \
    --variant "$gvcf_list_file" \
    -O "$combined_gvcf" \
    --tmp-dir "$output_dir/tmp" \
    --read-index 0 >> "$log_file" 2>&1

if [ $? -ne 0 ]; then
    echo "CombineGVCFs 失败" | tee -a "$log_file"
    exit 1
fi

echo "CombineGVCFs 成功完成" >> "$log_file"

# 创建GVCF索引
echo "创建GVCF索引..." >> "$log_file"
$gatk_path --java-options "-Xmx${java_memory} -Xms80g" IndexFeatureFile \
    -F "$combined_gvcf" >> "$log_file" 2>&1

if [ $? -ne 0 ]; then
    echo "索引创建失败" | tee -a "$log_file"
    exit 1
fi

echo "GVCF索引创建完成" >> "$log_file"

# 基本验证
echo "验证合并后的GVCF文件..." >> "$log_file"
$gatk_path --java-options "-Xmx20g" ValidateVariants \
    -V "$combined_gvcf" \
    -R "$reference_genome" >> "$log_file" 2>&1

if [ $? -eq 0 ]; then
    echo "GVCF文件验证通过" >> "$log_file"
else
    echo "警告: GVCF文件验证失败，但文件可能仍可用" >> "$log_file"
fi

# 简单统计
echo "收集基本统计信息..." >> "$log_file"
variant_count=$($gatk_path --java-options "-Xmx${java_memory} -Xms80g" CountVariants \
    -V "$combined_gvcf" 2>/dev/null | head -1)

echo "合并GVCF统计:" >> "$log_file"
echo "样本数量: $gvcf_count" >> "$log_file"
echo "总变异位点数: ${variant_count:-未知}" >> "$log_file"
echo "输出文件: $combined_gvcf" >> "$log_file"

# 清理
rm -rf "$output_dir/tmp"

echo "GATK CombineGVCFs 完成时间: $(date)" >> "$log_file"

# 显示结果
echo "GATK CombineGVCFs处理完成!"
echo "成功合并 $gvcf_count 个GVCF文件"
echo "合并文件: $combined_gvcf"
echo "索引文件: ${combined_gvcf}.tbi"
echo "日志文件: $log_file"
