#!/bin/bash

# 设置目录和文件
input_bam_dir="/home/data/3.bwa_sam_bam/sorted_bam"
output_dir="/home/data/4.picard_dedup"
picard_jar_path="/home/vensin/software/picard.jar"

# 创建输出目录
mkdir -p "$output_dir/markdup" "$output_dir/markdup/picard_metrics" "$output_dir/tmp" "$output_dir/logs"

# 设置日志文件
log_file="$output_dir/picard_processing.log"
echo "Picard MarkDuplicates processing started at $(date)" > "$log_file"

# 检查输入BAM文件数量
echo "=== 检查输入BAM文件 ===" >> "$log_file"
bam_count=$(find "$input_bam_dir" -name '*.sorted.bam' | wc -l)
echo "找到 $bam_count 个BAM文件" >> "$log_file"

if [ "$bam_count" -eq 0 ]; then
    echo "错误: 未找到排序后的BAM文件" | tee -a "$log_file"
    exit 1
fi

# 显示前几个样本
echo "BAM文件示例:" >> "$log_file"
find "$input_bam_dir" -name '*.sorted.bam' | head -5 | xargs -n 1 basename >> "$log_file"

# 创建Picard处理脚本
PICARD_SCRIPT="$output_dir/picard_process_sample.sh"

cat > "$PICARD_SCRIPT" << EOF
#!/bin/bash
# Picard样本处理脚本

sorted_bam_path=\$1
output_dir=\$2
picard_jar_path=\$3

base_name=\$(basename "\$sorted_bam_path" .sorted.bam)

# 样本特定的日志
sample_log="$output_dir/logs/\${base_name}.log"

echo "开始Picard处理样本: \$base_name" > "\$sample_log"
echo "  输入BAM: \$(basename "\$sorted_bam_path")" >> "\$sample_log"

# 创建样本特定的临时目录
sample_tmp_dir="$output_dir/tmp/\$base_name"
mkdir -p "\$sample_tmp_dir"

# 执行 Picard MarkDuplicates
echo "  运行MarkDuplicates..." >> "\$sample_log"
java -Xmx20g -jar "\$picard_jar_path" MarkDuplicates \\
    -I "\$sorted_bam_path" \\
    -O "$output_dir/markdup/\${base_name}.markdup.bam" \\
    -M "$output_dir/markdup/picard_metrics/\${base_name}.metrics.txt" \\
    --REMOVE_DUPLICATES false \\
    --ASSUME_SORTED true \\
    --VALIDATION_STRINGENCY LENIENT \\
    --TMP_DIR "\$sample_tmp_dir" 2>> "\$sample_log"

if [ \$? -eq 0 ]; then
    echo "  MarkDuplicates 成功: \$base_name" >> "\$sample_log"
    
    # 为标记重复后的BAM文件创建索引
    echo "  创建BAM索引..." >> "\$sample_log"
    samtools index "$output_dir/markdup/\${base_name}.markdup.bam" 2>> "\$sample_log"
    
    if [ \$? -eq 0 ]; then
        echo "  BAM索引创建成功: \$base_name" >> "\$sample_log"
        echo "\$base_name: 成功" >> "$output_dir/picard_processing.log"
    else
        echo "  BAM索引创建失败: \$base_name" >> "\$sample_log"
        echo "\$base_name: 索引失败" >> "$output_dir/picard_processing.log"
    fi
    
    # 清理临时文件
    rm -rf "\$sample_tmp_dir"
else
    echo "  MarkDuplicates 失败: \$base_name" >> "\$sample_log"
    echo "\$base_name: 失败" >> "$output_dir/picard_processing.log"
    # 清理临时文件
    rm -rf "\$sample_tmp_dir"
    exit 1
fi

echo "  完成: \$base_name" >> "\$sample_log"
EOF

chmod +x "$PICARD_SCRIPT"

# 处理所有样本 - 使用 GNU Parallel 进行并行处理
echo "=== 开始Picard MarkDuplicates处理 ===" >> "$log_file"

# 检查是否安装了 GNU Parallel
if command -v parallel >/dev/null 2>&1; then
    echo "使用 GNU Parallel 并行处理任务 ..." >> "$log_file"
    
    # 创建任务列表文件
    task_file="$output_dir/task_list.txt"
    find "$input_bam_dir" -name '*.sorted.bam' | sort > "$task_file"
    
    # 使用 parallel 并行处理
    cat "$task_file" | parallel -j 5 --joblog "$output_dir/parallel_jobs.log" \
        "$PICARD_SCRIPT {} $output_dir $picard_jar_path"
    
    # 计算成功数量
    successful_jobs=$(grep -c ": 成功" "$log_file" 2>/dev/null || echo 0)
else
    echo "GNU Parallel 未安装，使用串行处理..." >> "$log_file"
    
    # 串行处理
    success_count=0
    while IFS= read -r sorted_bam; do
        if [ -n "$sorted_bam" ]; then
            echo "处理样本: $(basename "$sorted_bam" .sorted.bam)" >> "$log_file"
            if "$PICARD_SCRIPT" "$sorted_bam" "$output_dir" "$picard_jar_path"; then
                ((success_count++))
            fi
        fi
    done < <(find "$input_bam_dir" -name '*.sorted.bam' | sort)
    
    successful_jobs=$success_count
fi

# 清理临时脚本
rm -f "$PICARD_SCRIPT"

# 修复算术运算 - 确保变量是数字
successful_jobs=${successful_jobs:-0}
bam_count=${bam_count:-0}
failed_count=$((bam_count - successful_jobs))

# 检查处理结果
echo "处理完成: $successful_jobs/$bam_count 个样本成功处理" >> "$log_file"

# 显示处理摘要
echo "=== Picard处理摘要 ===" >> "$log_file"
echo "成功: $successful_jobs" >> "$log_file"
echo "失败: $failed_count" >> "$log_file"

# 生成重复标记统计汇总
echo "=== 重复标记统计汇总 ===" >> "$log_file"
for metrics_file in "$output_dir/markdup/picard_metrics"/*.metrics.txt; do
    if [ -f "$metrics_file" ]; then
        sample_name=$(basename "$metrics_file" .metrics.txt)
        echo "样本: $sample_name" >> "$log_file"
        
        # 提取关键统计信息
        if grep -q "LIBRARY" "$metrics_file"; then
            echo "  重复统计:" >> "$log_file"
            grep -A 2 "LIBRARY" "$metrics_file" | tail -1 >> "$log_file"
            
            # 提取关键指标
            unpaired_reads=$(grep -A 2 "LIBRARY" "$metrics_file" | tail -1 | awk '{print $3}')
            read_pairs=$(grep -A 2 "LIBRARY" "$metrics_file" | tail -1 | awk '{print $4}')
            unpaired_duplicates=$(grep -A 2 "LIBRARY" "$metrics_file" | tail -1 | awk '{print $5}')
            paired_duplicates=$(grep -A 2 "LIBRARY" "$metrics_file" | tail -1 | awk '{print $6}')
            paired_optical_duplicates=$(grep -A 2 "LIBRARY" "$metrics_file" | tail -1 | awk '{print $7}')
            
            total_reads=$((unpaired_reads + read_pairs * 2))
            duplicate_reads=$((unpaired_duplicates + paired_duplicates * 2))
            duplication_rate=$(echo "scale=4; $duplicate_reads / $total_reads" | bc 2>/dev/null || echo "N/A")
            
            echo "  总读段数: $total_reads" >> "$log_file"
            echo "  重复读段数: $duplicate_reads" >> "$log_file"
            echo "  重复率: $duplication_rate" >> "$log_file"
        fi
        echo "---" >> "$log_file"
    fi
done

echo "Picard MarkDuplicates processing completed at $(date)" >> "$log_file"

# 在终端显示最终结果
echo "Picard MarkDuplicates处理完成!"
echo "成功处理: $successful_jobs/$bam_count 个样本"
echo "详细日志: $log_file"
echo "标记重复BAM文件: $output_dir/markdup/"
echo "统计文件: $output_dir/markdup/picard_metrics/"
