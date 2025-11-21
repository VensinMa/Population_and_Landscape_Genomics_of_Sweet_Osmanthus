#!/bin/bash

# 设置目录和文件
input_bam_dir="/home/data/4.picard_dedup/markdup"  # Picard去重后的BAM文件
output_dir="/home/data/5.gatk_haplotypecaller"
reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
gatk_path="/home/vensin/software/gatk-4.6.2.0/gatk" 

# 创建输出目录
mkdir -p "$output_dir/gvcf" "$output_dir/tmp" "$output_dir/logs" "$output_dir/metrics"

# 设置日志文件
log_file="$output_dir/haplotypecaller_processing.log"
echo "GATK HaplotypeCaller processing started at $(date)" > "$log_file"

# 检查输入BAM文件数量
echo "=== 检查输入BAM文件 ===" >> "$log_file"
bam_count=$(find "$input_bam_dir" -name '*.markdup.bam' | wc -l)
echo "找到 $bam_count 个标记重复的BAM文件" >> "$log_file"

if [ "$bam_count" -eq 0 ]; then
    echo "错误: 未找到标记重复的BAM文件" | tee -a "$log_file"
    exit 1
fi

# 显示前几个样本
echo "BAM文件示例:" >> "$log_file"
find "$input_bam_dir" -name '*.markdup.bam' | head -5 | xargs -n 1 basename >> "$log_file"

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

# 创建HaplotypeCaller处理脚本
HC_SCRIPT="$output_dir/haplotypecaller_process_sample.sh"

cat > "$HC_SCRIPT" << 'EOF'
#!/bin/bash
# HaplotypeCaller样本处理脚本

input_bam=$1
output_dir=$2
reference_genome=$3
gatk_path=$4

base_name=$(basename "$input_bam" .markdup.bam)

# 样本特定的日志
sample_log="$output_dir/logs/${base_name}.log"

echo "开始HaplotypeCaller处理样本: $base_name" > "$sample_log"
echo "  输入BAM: $(basename "$input_bam")" >> "$sample_log"

# 创建样本特定的临时目录
sample_tmp_dir="$output_dir/tmp/$base_name"
mkdir -p "$sample_tmp_dir"

# 输出文件
gvcf_output="$output_dir/gvcf/${base_name}.g.vcf.gz"

# 执行 GATK HaplotypeCaller
echo "  运行HaplotypeCaller..." >> "$sample_log"
$gatk_path --java-options "-Xmx20g -Xms4g" HaplotypeCaller \
    -R "$reference_genome" \
    -I "$input_bam" \
    -O "$gvcf_output" \
    -ERC GVCF \
    --tmp-dir "$sample_tmp_dir" \
    --native-pair-hmm-threads 4 \
    --sample-ploidy 2 >> "$sample_log"

hc_exit_code=$?

if [ $hc_exit_code -eq 0 ]; then
    echo "  HaplotypeCaller 成功: $base_name" >> "$sample_log"
    
    # 验证输出的GVCF文件
    if [ -f "$gvcf_output" ]; then
        echo "  GVCF文件生成成功: $(basename "$gvcf_output")" >> "$sample_log"
        echo "$base_name: 成功" >> "$output_dir/haplotypecaller_processing.log"
        
        # 可选：收集一些基本统计
        echo "  收集GVCF统计..." >> "$sample_log"
        $gatk_path --java-options "-Xmx4g" CountVariants \
            -V "$gvcf_output" 2>> "$sample_log" | head -5 >> "$sample_log"
    else
        echo "  GVCF文件未生成: $base_name" >> "$sample_log"
        echo "$base_name: GVCF文件缺失" >> "$output_dir/haplotypecaller_processing.log"
        hc_exit_code=1
    fi
else
    echo "  HaplotypeCaller 失败 (退出码: $hc_exit_code): $base_name" >> "$sample_log"
    echo "$base_name: 失败" >> "$output_dir/haplotypecaller_processing.log"
fi

# 清理临时文件
rm -rf "$sample_tmp_dir"

echo "  处理完成: $base_name" >> "$sample_log"
exit $hc_exit_code
EOF

chmod +x "$HC_SCRIPT"

# 处理所有样本 - 使用 GNU Parallel 进行并行处理
echo "=== 开始GATK HaplotypeCaller处理 ===" >> "$log_file"

# 检查是否安装了 GNU Parallel
if command -v parallel >/dev/null 2>&1; then
    echo "使用 GNU Parallel 并行处理 (2个并行任务，因为HaplotypeCaller比较耗资源)..." >> "$log_file"
    
    # 创建任务列表文件
    task_file="$output_dir/task_list.txt"
    find "$input_bam_dir" -name '*.markdup.bam' | sort > "$task_file"
    
    # 使用 parallel 并行处理
    cat "$task_file" | parallel -j 2 --joblog "$output_dir/parallel_jobs.log" \
        "$HC_SCRIPT {} $output_dir $reference_genome $gatk_path"
    
    # 计算成功数量
    successful_jobs=$(grep -c ": 成功" "$log_file" 2>/dev/null || echo 0)
else
    echo "GNU Parallel 未安装，使用串行处理..." >> "$log_file"
    
    # 串行处理
    success_count=0
    processed_count=0
    while IFS= read -r input_bam; do
        if [ -n "$input_bam" ]; then
            ((processed_count++))
            echo "处理样本 $processed_count/$bam_count: $(basename "$input_bam" .markdup.bam)" | tee -a "$log_file"
            if "$HC_SCRIPT" "$input_bam" "$output_dir" "$reference_genome" "$gatk_path"; then
                ((success_count++))
            fi
        fi
    done < <(find "$input_bam_dir" -name '*.markdup.bam' | sort)
    
    successful_jobs=$success_count
fi

# 清理临时脚本
rm -f "$HC_SCRIPT"

# 统计结果
successful_jobs=${successful_jobs:-0}
bam_count=${bam_count:-0}
failed_count=$((bam_count - successful_jobs))

echo "处理完成: $successful_jobs/$bam_count 个样本成功处理" >> "$log_file"

# 显示处理摘要
echo "=== HaplotypeCaller处理摘要 ===" >> "$log_file"
echo "成功: $successful_jobs" >> "$log_file"
echo "失败: $failed_count" >> "$log_file"

# 生成GVCF文件统计汇总
echo "=== GVCF文件统计汇总 ===" >> "$log_file"
for gvcf_file in "$output_dir/gvcf"/*.g.vcf.gz; do
    if [ -f "$gvcf_file" ]; then
        sample_name=$(basename "$gvcf_file" .g.vcf.gz)
        echo "样本: $sample_name" >> "$log_file"
        
        # 检查文件完整性
        if $gatk_path ValidateVariants -V "$gvcf_file" 2>/dev/null; then
            echo "  文件验证: 通过" >> "$log_file"
            
            # 获取变体数量统计
            variant_count=$($gatk_path --java-options "-Xmx4g" CountVariants -V "$gvcf_file" 2>/dev/null | awk '/CountVariants/{getline; print $1}')
            if [ -n "$variant_count" ]; then
                echo "  变体数量: $variant_count" >> "$log_file"
            fi
        else
            echo "  文件验证: 失败" >> "$log_file"
        fi
        echo "---" >> "$log_file"
    fi
done

echo "GATK HaplotypeCaller processing completed at $(date)" >> "$log_file"

# 在终端显示最终结果
echo "GATK HaplotypeCaller处理完成!"
echo "成功处理: $successful_jobs/$bam_count 个样本"
echo "详细日志: $log_file"
echo "GVCF文件: $output_dir/gvcf/"
echo "各样本日志: $output_dir/logs/"
