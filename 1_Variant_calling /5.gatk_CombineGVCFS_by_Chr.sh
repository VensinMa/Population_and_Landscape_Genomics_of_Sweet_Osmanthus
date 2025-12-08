#!/bin/bash

# 设置目录和文件
input_gvcf_dir="/home/vensin/workspace/snpcalling_wild/5.gatk_haplotypecaller/gvcf"
output_dir="/home/data/6.gatk_combinegvcf_per_chr"
reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
gatk_path="/home/vensin/software/gatk-4.6.2.0/gatk"

# 默认并行任务数
PARALLEL_JOBS=6

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -j|--jobs)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        -h|--help)
            echo "用法: $0 [-j 并行任务数]"
            echo "   -j, --jobs     指定并行任务数（默认: 6）"
            exit 0
            ;;
        *)
            echo "未知参数: $1"
            echo "使用 -h 查看帮助"
            exit 1
            ;;
    esac
done

# 验证并行任务数
if ! [[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || [ "$PARALLEL_JOBS" -lt 1 ]; then
    echo "错误: 并行任务数必须是正整数"
    exit 1
fi

echo "设置并行任务数为: $PARALLEL_JOBS"

# 创建输出目录
mkdir -p "$output_dir" "$output_dir/per_chromosome" "$output_dir/tmp" "$output_dir/logs"

# 设置日志文件
log_file="$output_dir/combinegvcfs.log"
echo "GATK CombineGVCFs 按染色体合并开始时间: $(date)" > "$log_file"

# 检查输入GVCF文件
gvcf_files=($(find "$input_gvcf_dir" -name '*.g.vcf.gz'))
gvcf_count=${#gvcf_files[@]}

if [ "$gvcf_count" -eq 0 ]; then
    echo "错误: 未找到GVCF文件" | tee -a "$log_file"
    exit 1
fi

echo "找到 $gvcf_count 个GVCF文件" >> "$log_file"

# 显示前几个样本
echo "GVCF文件示例:" >> "$log_file"
find "$input_gvcf_dir" -name '*.g.vcf.gz' | head -1000 | xargs -n 1 basename >> "$log_file"

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
gvcf_list_file="$output_dir/gvcf_files.list"
printf "%s\n" "${gvcf_files[@]}" | sort > "$gvcf_list_file"

echo "GVCF文件列表已创建: $gvcf_list_file" >> "$log_file"
echo "包含以下文件:" >> "$log_file"
head -5 "$gvcf_list_file" >> "$log_file"
echo "[...]" >> "$log_file"

# 染色体列表
chromosomes=("Chr01" "Chr02" "Chr03" "Chr04" "Chr05" "Chr06" "Chr07" "Chr08" "Chr09" "Chr10"
             "Chr11" "Chr12" "Chr13" "Chr14" "Chr15" "Chr16" "Chr17" "Chr18" "Chr19" "Chr20"
             "Chr21" "Chr22" "Chr23")

echo "开始按染色体合并GVCF文件..." >> "$log_file"
echo "共 ${#chromosomes[@]} 条染色体需要处理" >> "$log_file"

# 创建CombineGVCFs染色体处理脚本（参考成功脚本的方法）
COMBINE_SCRIPT="$output_dir/combinegvcfs_process_chromosome.sh"

cat > "$COMBINE_SCRIPT" << 'EOF'
#!/bin/bash
# CombineGVCFs染色体处理脚本

chrom=$1
output_dir=$2
reference_genome=$3
gatk_path=$4
gvcf_list_file=$5

# 染色体特定的日志
chrom_log="$output_dir/logs/combine_${chrom}.log"

echo "开始CombineGVCFs处理染色体: $chrom" > "$chrom_log"
echo "  染色体: $chrom" >> "$chrom_log"
echo "  输入GVCF文件列表: $(basename "$gvcf_list_file")" >> "$chrom_log"

# 创建染色体特定的临时目录
chrom_tmp_dir="$output_dir/tmp/$chrom"
mkdir -p "$chrom_tmp_dir"

# 输出文件
gvcf_output="$output_dir/per_chromosome/combined.${chrom}.g.vcf.gz"

# 执行 GATK CombineGVCFs（按染色体）
echo "  运行CombineGVCFs (染色体: $chrom)..." >> "$chrom_log"
$gatk_path --java-options "-Xmx40g -Xms10g" CombineGVCFs \
    -R "$reference_genome" \
    --variant "$gvcf_list_file" \
    -L "$chrom" \
    -O "$gvcf_output" \
    --tmp-dir "$chrom_tmp_dir" \
    --read-index 0 >> "$chrom_log" 2>&1

combine_exit_code=$?

if [ $combine_exit_code -eq 0 ]; then
    echo "  CombineGVCFs 成功: 染色体 $chrom" >> "$chrom_log"
    
    # 验证输出的GVCF文件
    if [ -f "$gvcf_output" ]; then
        echo "  GVCF文件生成成功: $(basename "$gvcf_output")" >> "$chrom_log"
        echo "$chrom: 成功" >> "$output_dir/combinegvcfs.log"
        
        # 创建GVCF索引
        echo "  创建GVCF索引..." >> "$chrom_log"
        $gatk_path --java-options "-Xmx40g -Xms10g" IndexFeatureFile \
            -F "$gvcf_output" >> "$chrom_log" 2>&1
            
        if [ $? -eq 0 ]; then
            echo "  索引创建成功: 染色体 $chrom" >> "$chrom_log"
        else
            echo "  索引创建失败: 染色体 $chrom" >> "$chrom_log"
            combine_exit_code=1
        fi
    else
        echo "  GVCF文件未生成: 染色体 $chrom" >> "$chrom_log"
        echo "$chrom: GVCF文件缺失" >> "$output_dir/combinegvcfs.log"
        combine_exit_code=1
    fi
else
    echo "  CombineGVCFs 失败 (退出码: $combine_exit_code): 染色体 $chrom" >> "$chrom_log"
    echo "$chrom: 失败" >> "$output_dir/combinegvcfs.log"
fi

# 清理临时文件
rm -rf "$chrom_tmp_dir"

echo "  处理完成: 染色体 $chrom" >> "$chrom_log"
exit $combine_exit_code
EOF

chmod +x "$COMBINE_SCRIPT"

# 处理所有染色体 - 使用 GNU Parallel 进行并行处理
echo "=== 开始GATK CombineGVCFs按染色体处理 ===" >> "$log_file"

# 检查是否安装了 GNU Parallel
if command -v parallel >/dev/null 2>&1; then
    echo "使用 GNU Parallel 并行处理 (${PARALLEL_JOBS}个并行任务)..." >> "$log_file"
    
    # 创建任务列表文件（染色体列表）
    task_file="$output_dir/chromosome_list.txt"
    printf "%s\n" "${chromosomes[@]}" > "$task_file"
    
    echo "染色体列表:" >> "$log_file"
    cat "$task_file" >> "$log_file"
    
    # 使用 parallel 并行处理
    cat "$task_file" | parallel -j $PARALLEL_JOBS --joblog "$output_dir/parallel_jobs.log" \
        "$COMBINE_SCRIPT {} $output_dir $reference_genome $gatk_path $gvcf_list_file"
    
    # 计算成功数量
    successful_jobs=$(grep -c ": 成功" "$output_dir/combinegvcfs.log" 2>/dev/null || echo 0)
else
    echo "GNU Parallel 未安装，使用串行处理..." >> "$log_file"
    
    # 串行处理
    success_count=0
    processed_count=0
    for chrom in "${chromosomes[@]}"; do
        ((processed_count++))
        echo "处理染色体 $processed_count/${#chromosomes[@]}: $chrom" | tee -a "$log_file"
        if "$COMBINE_SCRIPT" "$chrom" "$output_dir" "$reference_genome" "$gatk_path" "$gvcf_list_file"; then
            ((success_count++))
        fi
    done
    
    successful_jobs=$success_count
fi

# 清理临时脚本
rm -f "$COMBINE_SCRIPT"

# 统计结果
successful_jobs=${successful_jobs:-0}
chrom_count=${#chromosomes[@]}
failed_count=$((chrom_count - successful_jobs))

echo "处理完成: $successful_jobs/$chrom_count 条染色体成功处理" >> "$log_file"

# 显示处理摘要
echo "=== CombineGVCFs按染色体处理摘要 ===" >> "$log_file"
echo "成功: $successful_jobs" >> "$log_file"
echo "失败: $failed_count" >> "$log_file"

# 生成染色体文件统计汇总
echo "=== 染色体GVCF文件统计汇总 ===" >> "$log_file"
for chrom in "${chromosomes[@]}"; do
    gvcf_file="$output_dir/per_chromosome/combined.${chrom}.g.vcf.gz"
    if [ -f "$gvcf_file" ] && [ -f "${gvcf_file}.tbi" ]; then
        echo "染色体: $chrom" >> "$log_file"
        
        # 获取文件大小
        file_size=$(ls -lh "$gvcf_file" | awk '{print $5}')
        echo "  文件大小: $file_size" >> "$log_file"
        
        # 可选：获取变体数量统计
        echo "  文件位置: $gvcf_file" >> "$log_file"
    else
        echo "染色体: $chrom - 文件缺失或索引未创建" >> "$log_file"
    fi
done

# 生成染色体文件列表（用于后续分析）
chrom_list_file="$output_dir/chromosome_files.list"
> "$chrom_list_file"
for chrom in "${chromosomes[@]}"; do
    chrom_file="$output_dir/per_chromosome/combined.${chrom}.g.vcf.gz"
    if [ -f "$chrom_file" ] && [ -f "${chrom_file}.tbi" ]; then
        echo "$chrom_file" >> "$chrom_list_file"
    fi
done

echo "染色体文件列表已创建: $chrom_list_file" >> "$log_file"
echo "包含 $(wc -l < "$chrom_list_file") 个染色体文件" >> "$log_file"

# 如果需要，可以合并所有染色体的文件
if [ $successful_jobs -eq $chrom_count ]; then
    echo "所有染色体合并成功，可以创建最终合并文件（可选）..." >> "$log_file"
    
    final_combined="$output_dir/combined.all_chromosomes.g.vcf.gz"
    
    echo "合并所有染色体文件为单个文件..." >> "$log_file"
    $gatk_path --java-options "-Xmx100g" MergeVcfs \
        -R "$reference_genome" \
        --variant "$chrom_list_file" \
        -O "$final_combined" \
        --tmp-dir "$output_dir/tmp/final_merge" >> "$log_file" 2>&1
        
    if [ $? -eq 0 ]; then
        echo "最终合并文件创建成功: $final_combined" >> "$log_file"
        
        # 创建最终索引
        $gatk_path --java-options "-Xmx100g" IndexFeatureFile \
            -F "$final_combined" >> "$log_file" 2>&1
    else
        echo "最终合并文件创建失败" >> "$log_file"
    fi
else
    echo "有 $failed_count 条染色体合并失败，跳过最终合并步骤" >> "$log_file"
fi

echo "GATK CombineGVCFs processing completed at $(date)" >> "$log_file"

# 在终端显示最终结果
echo ""
echo "GATK CombineGVCFs 按染色体合并完成!"
echo "成功合并: $successful_jobs/$chrom_count 条染色体"
echo "详细日志: $log_file"
echo "染色体GVCF文件: $output_dir/per_chromosome/"
echo "染色体文件列表: $chrom_list_file"
echo "各染色体日志: $output_dir/logs/"

if [ $failed_count -gt 0 ]; then
    echo "注意: 有 $failed_count 条染色体合并失败，请检查日志文件"
    echo "失败的染色体日志: $output_dir/logs/combine_*.log"
fi
