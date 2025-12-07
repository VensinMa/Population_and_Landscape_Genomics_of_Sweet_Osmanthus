#!/bin/bash

# 设置目录和文件
input_gvcf_dir="/home/vensin/workspace/snpcalling_wild/5.gatk_haplotypecaller/gvcf"  # HaplotypeCaller生成的GVCF文件
output_dir="/home/data/6.gatk_combinegvcfs"
reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
gatk_path="/home/vensin/software/gatk-4.6.2.0/gatk"

# 创建输出目录
mkdir -p "$output_dir/combined_gvcf" "$output_dir/tmp" "$output_dir/logs" "$output_dir/gvcf_lists"

# 设置日志文件
log_file="$output_dir/combinegvcfs_processing.log"
echo "GATK CombineGVCFs processing started at $(date)" > "$log_file"

# 检查输入GVCF文件数量
echo "=== 检查输入GVCF文件 ===" >> "$log_file"
gvcf_count=$(find "$input_gvcf_dir" -name '*.g.vcf.gz' | wc -l)
echo "找到 $gvcf_count 个GVCF文件" >> "$log_file"

if [ "$gvcf_count" -eq 0 ]; then
    echo "错误: 未找到GVCF文件" | tee -a "$log_file"
    exit 1
fi

# 显示前几个GVCF文件
echo "GVCF文件示例:" >> "$log_file"
find "$input_gvcf_dir" -name '*.g.vcf.gz' | head -5 | xargs -n 1 basename >> "$log_file"

# 检查参考基因组和相关文件
echo "=== 检查参考基因组和相关文件 ===" >> "$log_file"
if [ ! -f "$reference_genome" ]; then
    echo "错误: 参考基因组文件不存在: $reference_genome" | tee -a "$log_file"
    exit 1
fi

if [ ! -f "$reference_genome.fai" ]; then
    echo "创建FASTA索引..." >> "$log_file"
    samtools faidx "$reference_genome" >> "$log_file" 2>&1
fi

dict_file="${reference_genome%.*}.dict"
if [ ! -f "$dict_file" ]; then
    echo "创建序列字典..." >> "$log_file"
    $gatk_path CreateSequenceDictionary -R "$reference_genome" -O "$dict_file" >> "$log_file" 2>&1
fi

# 创建GVCF文件列表
echo "=== 创建GVCF文件列表 ===" >> "$log_file"
gvcf_list_file="$output_dir/gvcf_lists/gvcf_files.list"
find "$input_gvcf_dir" -name '*.g.vcf.gz' | sort > "$gvcf_list_file"

echo "GVCF文件列表已创建: $gvcf_list_file" >> "$log_file"
echo "包含以下文件:" >> "$log_file"
cat "$gvcf_list_file" >> "$log_file"

# 合并GVCF文件的选项
# 根据样本数量调整内存
if [ "$gvcf_count" -gt 200 ]; then
    java_memory="100g"
elif [ "$gvcf_count" -gt 100 ]; then
    java_memory="100g"
else
    java_memory="100g"
fi

echo "根据 $gvcf_count 个样本，设置内存为: $java_memory" >> "$log_file"

# 合并GVCF文件
echo "=== 开始合并GVCF文件 ===" >> "$log_file"
combined_gvcf="$output_dir/combined_gvcf/combined.g.vcf.gz"

# 创建临时目录
combined_tmp_dir="$output_dir/tmp/combinegvcfs"
mkdir -p "$combined_tmp_dir"

# 执行 GATK CombineGVCFs
echo "运行CombineGVCFs..." >> "$log_file"
$gatk_path --java-options "-Xmx${java_memory} -Xms80g" CombineGVCFs \
    -R "$reference_genome" \
    --variant "$gvcf_list_file" \
    -O "$combined_gvcf" \
    --tmp-dir "$combined_tmp_dir" \
    --read-index 0 2>> "$log_file"

combine_exit_code=$?

if [ $combine_exit_code -eq 0 ]; then
    echo "CombineGVCFs 成功" >> "$log_file"
    
    # 验证合并后的GVCF文件
    echo "验证合并后的GVCF文件..." >> "$log_file"
    $gatk_path --java-options "-Xmx100g -Xms20g" ValidateVariants \
        -V "$combined_gvcf" \
        -R "$reference_genome" 2>> "$log_file"
    
    if [ $? -eq 0 ]; then
        echo "GVCF文件验证通过" >> "$log_file"
    else
        echo "GVCF文件验证失败" >> "$log_file"
    fi
    
    # 收集合并后的GVCF统计信息
    echo "收集合并GVCF统计信息..." >> "$log_file"
    stats_file="$output_dir/combined_gvcf/combined_gvcf_stats.txt"
    $gatk_path --java-options "-Xmx100g -Xms20g" CountVariants \
        -V "$combined_gvcf" > "$stats_file" 2>> "$log_file"
    
    # 创建合并GVCF的索引（如果未自动创建）
    echo "创建合并GVCF索引..." >> "$log_file"
    if [ ! -f "${combined_gvcf}.tbi" ]; then
        $gatk_path --java-options "-Xmx100g -Xms20g" IndexFeatureFile \
            -F "$combined_gvcf" 2>> "$log_file"
    fi
    
    # 显示基本统计信息
    echo "=== 合并GVCF统计摘要 ===" >> "$log_file"
    echo "样本数量: $gvcf_count" >> "$log_file"
    
    if [ -f "$stats_file" ]; then
        variant_count=$(grep -E "^[0-9]+$" "$stats_file" | head -1)
        if [ -n "$variant_count" ]; then
            echo "总变异位点数: $variant_count" >> "$log_file"
        fi
    fi
    
    # 获取文件大小信息
    combined_size=$(ls -lh "$combined_gvcf" | awk '{print $5}')
    echo "合并GVCF文件大小: $combined_size" >> "$log_file"
    
    # 检查每个染色体的变异数量（可选）
    echo "按染色体统计变异数量:" >> "$log_file"
    $gatk_path --java-options "-Xmx100g -Xms20g" CountVariants \
        -V "$combined_gvcf" \
        -L "Chr01" -L "Chr02" -L "Chr03" -L "Chr04" -L "Chr05" -L "Chr06" -L "Chr07" -L "Chr08" -L "Chr09" -L "Chr10" \
        -L "Chr11" -L "Chr12" -L "Chr13" -L "Chr14" -L "Chr15" -L "Chr16" -L "Chr17" -L "Chr18" -L "Chr19" -L "Chr20" \
        -L "Chr21" -L "Chr22" -L "Chr23"  2>/dev/null | \
        grep -E "^Chr" >> "$log_file" 2>/dev/null || echo "  无法按染色体统计" >> "$log_file"
    
    echo "$gvcf_count 个样本合并成功" >> "$output_dir/combinegvcfs_processing.log"
else
    echo "CombineGVCFs 失败 (退出码: $combine_exit_code)" >> "$log_file"
    echo "$gvcf_count 个样本合并失败" >> "$output_dir/combinegvcfs_processing.log"
fi

# 清理临时文件
rm -rf "$combined_tmp_dir"

# 统计结果
if [ $combine_exit_code -eq 0 ]; then
    successful_merge=1
else
    successful_merge=0
fi

echo "处理完成: $successful_merge/1 个合并任务" >> "$log_file"

# 显示处理摘要
echo "=== CombineGVCFs处理摘要 ===" >> "$log_file"
if [ $successful_merge -eq 1 ]; then
    echo "状态: 成功" >> "$log_file"
    echo "合并后的GVCF文件: $combined_gvcf" >> "$log_file"
else
    echo "状态: 失败" >> "$log_file"
fi

echo "输入GVCF文件数: $gvcf_count" >> "$log_file"
echo "输出目录: $output_dir" >> "$log_file"

echo "GATK CombineGVCFs processing completed at $(date)" >> "$log_file"

# 在终端显示最终结果
if [ $successful_merge -eq 1 ]; then
    echo "GATK CombineGVCFs处理完成!"
    echo "成功合并: $gvcf_count 个样本的GVCF文件"
    echo "合并后的GVCF文件: $combined_gvcf"
    if [ -f "$stats_file" ]; then
        variant_count=$(grep -E "^[0-9]+$" "$stats_file" | head -1)
        if [ -n "$variant_count" ]; then
            echo "总变异位点数: $variant_count"
        fi
    fi
else
    echo "GATK CombineGVCFs处理失败!"
    echo "请查看日志文件: $log_file"
fi
echo "详细日志: $log_file"
echo "GVCF文件列表: $gvcf_list_file"
