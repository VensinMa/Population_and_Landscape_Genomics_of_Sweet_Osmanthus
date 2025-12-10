#!/bin/bash

# 设置目录和文件
input_gvcf_dir="/home/data/6.gatk_combinegvcf_per_chr/per_chromosome"  # 上一步按染色体合并的GVCF文件
output_dir="/home/data/7.gatk_genotypegvcfs_per_chr"
reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
gatk_path="/home/vensin/software/gatk-4.6.2.0/gatk"

# 默认并行任务数
PARALLEL_JOBS=12

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -j|--jobs)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        -h|--help)
            echo "用法: $0 [-j 并行任务数]"
            echo "   -j, --jobs     指定并行任务数（默认: 12）"
            echo ""
            echo "输入文件目录: $input_gvcf_dir"
            echo "输出文件目录: $output_dir"
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
log_file="$output_dir/genotypegvcfs.log"
echo "GATK GenotypeGVCFs 按染色体基因分型开始时间: $(date)" > "$log_file"

# 检查输入GVCF文件（上一步按染色体合并的GVCF文件）
echo "=== 检查输入GVCF文件 ===" >> "$log_file"
gvcf_count=$(find "$input_gvcf_dir" -name 'combined.*.g.vcf.gz' | wc -l)

if [ "$gvcf_count" -eq 0 ]; then
    echo "错误: 未找到按染色体合并的GVCF文件" | tee -a "$log_file"
    echo "请检查输入目录: $input_gvcf_dir" | tee -a "$log_file"
    echo "预期文件格式: combined.Chr01.g.vcf.gz, combined.Chr02.g.vcf.gz 等" | tee -a "$log_file"
    exit 1
fi

echo "找到 $gvcf_count 个按染色体合并的GVCF文件" >> "$log_file"

# 显示前几个染色体文件
echo "染色体GVCF文件示例:" >> "$log_file"
find "$input_gvcf_dir" -name 'combined.*.g.vcf.gz' | head -10 | xargs -n 1 basename >> "$log_file"

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

# 获取染色体列表（从输入文件名中提取）
echo "=== 提取染色体列表 ===" >> "$log_file"
chromosomes=()
for gvcf_file in "$input_gvcf_dir"/combined.*.g.vcf.gz; do
    if [ -f "$gvcf_file" ]; then
        # 提取染色体名，例如从 combined.Chr01.g.vcf.gz 提取 Chr01
        chrom=$(basename "$gvcf_file" .g.vcf.gz | sed 's/combined\.//')
        chromosomes+=("$chrom")
    fi
done

# 如果没有找到染色体，使用默认列表
if [ ${#chromosomes[@]} -eq 0 ]; then
    echo "警告: 从文件名中未提取到染色体，使用默认染色体列表" >> "$log_file"
    chromosomes=("Chr01" "Chr02" "Chr03" "Chr04" "Chr05" "Chr06" "Chr07" "Chr08" "Chr09" "Chr10"
                 "Chr11" "Chr12" "Chr13" "Chr14" "Chr15" "Chr16" "Chr17" "Chr18" "Chr19" "Chr20"
                 "Chr21" "Chr22" "Chr23")
fi

# 排序染色体列表
IFS=$'\n' chromosomes=($(sort <<<"${chromosomes[*]}"))
unset IFS

echo "找到 ${#chromosomes[@]} 条染色体需要处理:" >> "$log_file"
printf "%s\n" "${chromosomes[@]}" >> "$log_file"

# 创建GenotypeGVCFs染色体处理脚本
GENOTYPE_SCRIPT="$output_dir/genotypegvcfs_process_chromosome.sh"

cat > "$GENOTYPE_SCRIPT" << 'EOF'
#!/bin/bash
# GenotypeGVCFs染色体处理脚本

chrom=$1
output_dir=$2
reference_genome=$3
gatk_path=$4
input_gvcf_dir=$5

# 染色体特定的日志
chrom_log="$output_dir/logs/genotype_${chrom}.log"

echo "开始GenotypeGVCFs处理染色体: $chrom" > "$chrom_log"
echo "  染色体: $chrom" >> "$chrom_log"

# 输入文件（上一步按染色体合并的GVCF文件）
input_gvcf="$input_gvcf_dir/combined.${chrom}.g.vcf.gz"

if [ ! -f "$input_gvcf" ]; then
    echo "错误: 输入GVCF文件不存在: $input_gvcf" >> "$chrom_log"
    echo "$chrom: 输入文件缺失" >> "$output_dir/genotypegvcfs.log"
    exit 1
fi

# 检查索引文件是否存在
if [ ! -f "${input_gvcf}.tbi" ]; then
    echo "警告: 输入GVCF索引文件不存在，尝试创建..." >> "$chrom_log"
    $gatk_path --java-options "-Xmx10g" IndexFeatureFile -I "$input_gvcf" >> "$chrom_log" 2>&1
    if [ $? -ne 0 ]; then
        echo "错误: 无法创建GVCF索引" >> "$chrom_log"
        echo "$chrom: 索引创建失败" >> "$output_dir/genotypegvcfs.log"
        exit 1
    fi
fi

# 创建染色体特定的临时目录
chrom_tmp_dir="$output_dir/tmp/$chrom"
mkdir -p "$chrom_tmp_dir"

# 输出文件
vcf_output="$output_dir/per_chromosome/genotyped.${chrom}.vcf.gz"

# 执行 GATK GenotypeGVCFs（按染色体）- 使用默认参数
echo "  运行GenotypeGVCFs (染色体: $chrom)..." >> "$chrom_log"
echo "  输入文件: $(basename "$input_gvcf")" >> "$chrom_log"
echo "  输出文件: $(basename "$vcf_output")" >> "$chrom_log"

$gatk_path --java-options "-Xmx40g -Xms10g" GenotypeGVCFs \
    -R "$reference_genome" \
    -V "$input_gvcf" \
    -O "$vcf_output" \
    --tmp-dir "$chrom_tmp_dir" >> "$chrom_log" 2>&1

genotype_exit_code=$?

if [ $genotype_exit_code -eq 0 ]; then
    echo "  GenotypeGVCFs 成功: 染色体 $chrom" >> "$chrom_log"
    
    # 验证输出的VCF文件
    if [ -f "$vcf_output" ]; then
        echo "  VCF文件生成成功: $(basename "$vcf_output")" >> "$chrom_log"
        echo "$chrom: 成功" >> "$output_dir/genotypegvcfs.log"
        
        # 创建VCF索引
        echo "  创建VCF索引..." >> "$chrom_log"
        $gatk_path --java-options "-Xmx10g" IndexFeatureFile \
            -I "$vcf_output" >> "$chrom_log" 2>&1
            
        if [ $? -eq 0 ]; then
            echo "  索引创建成功: 染色体 $chrom" >> "$chrom_log"
        else
            echo "  索引创建失败: 染色体 $chrom" >> "$chrom_log"
            genotype_exit_code=1
        fi
    else
        echo "  VCF文件未生成: 染色体 $chrom" >> "$chrom_log"
        echo "$chrom: VCF文件缺失" >> "$output_dir/genotypegvcfs.log"
        genotype_exit_code=1
    fi
else
    echo "  GenotypeGVCFs 失败 (退出码: $genotype_exit_code): 染色体 $chrom" >> "$chrom_log"
    echo "$chrom: 失败" >> "$output_dir/genotypegvcfs.log"
fi

# 清理临时文件
rm -rf "$chrom_tmp_dir"

echo "  处理完成: 染色体 $chrom" >> "$chrom_log"
exit $genotype_exit_code
EOF

chmod +x "$GENOTYPE_SCRIPT"

# 处理所有染色体 - 使用 GNU Parallel 进行并行处理
echo "=== 开始GATK GenotypeGVCFs按染色体基因分型 ===" >> "$log_file"

# 检查是否安装了 GNU Parallel
if command -v parallel >/dev/null 2>&1; then
    echo "使用 GNU Parallel 并行处理 (${PARALLEL_JOBS}个并行任务)..." >> "$log_file"
    
    # 创建任务列表文件（染色体列表）
    task_file="$output_dir/chromosome_list.txt"
    printf "%s\n" "${chromosomes[@]}" > "$task_file"
    
    echo "染色体列表:" >> "$log_file"
    cat "$task_file" >> "$log_file"
    
    echo "开始并行处理，请等待..." >> "$log_file"
    
    # 使用 parallel 并行处理
    cat "$task_file" | parallel -j $PARALLEL_JOBS --joblog "$output_dir/parallel_jobs.log" \
        --progress \
        --eta \
        "$GENOTYPE_SCRIPT {} $output_dir $reference_genome $gatk_path $input_gvcf_dir"
    
    # 计算成功数量
    successful_jobs=$(grep -c ": 成功" "$output_dir/genotypegvcfs.log" 2>/dev/null || echo 0)
    parallel_exit_code=$?
    
    if [ $parallel_exit_code -ne 0 ]; then
        echo "Parallel 执行过程中出现错误 (退出码: $parallel_exit_code)" >> "$log_file"
    fi
else
    echo "GNU Parallel 未安装，使用串行处理..." >> "$log_file"
    
    # 串行处理
    success_count=0
    processed_count=0
    for chrom in "${chromosomes[@]}"; do
        ((processed_count++))
        echo "处理染色体 $processed_count/${#chromosomes[@]}: $chrom" | tee -a "$log_file"
        if "$GENOTYPE_SCRIPT" "$chrom" "$output_dir" "$reference_genome" "$gatk_path" "$input_gvcf_dir"; then
            ((success_count++))
        fi
    done
    
    successful_jobs=$success_count
fi

# 清理临时脚本
rm -f "$GENOTYPE_SCRIPT"

# 统计结果
successful_jobs=${successful_jobs:-0}
chrom_count=${#chromosomes[@]}
failed_count=$((chrom_count - successful_jobs))

echo "" >> "$log_file"
echo "处理完成: $successful_jobs/$chrom_count 条染色体成功处理" >> "$log_file"

# 显示处理摘要
echo "=== GenotypeGVCFs按染色体处理摘要 ===" >> "$log_file"
echo "成功: $successful_jobs" >> "$log_file"
echo "失败: $failed_count" >> "$log_file"

if [ $failed_count -gt 0 ]; then
    echo "失败的染色体:" >> "$log_file"
    grep ": 失败" "$output_dir/genotypegvcfs.log" >> "$log_file"
fi

# 生成染色体文件统计汇总
echo "" >> "$log_file"
echo "=== 染色体VCF文件统计汇总 ===" >> "$log_file"
for chrom in "${chromosomes[@]}"; do
    vcf_file="$output_dir/per_chromosome/genotyped.${chrom}.vcf.gz"
    if [ -f "$vcf_file" ] && [ -f "${vcf_file}.tbi" ]; then
        echo "染色体: $chrom" >> "$log_file"
        
        # 获取文件大小
        file_size=$(ls -lh "$vcf_file" | awk '{print $5}')
        echo "  文件大小: $file_size" >> "$log_file"
        echo "  文件位置: $vcf_file" >> "$log_file"
        echo "---" >> "$log_file"
    else
        echo "染色体: $chrom - 文件缺失或索引未创建" >> "$log_file"
    fi
done

# 生成染色体文件列表（用于后续分析）
chrom_list_file="$output_dir/chromosome_files.list"
> "$chrom_list_file"
for chrom in "${chromosomes[@]}"; do
    chrom_file="$output_dir/per_chromosome/genotyped.${chrom}.vcf.gz"
    if [ -f "$chrom_file" ] && [ -f "${chrom_file}.tbi" ]; then
        echo "$chrom_file" >> "$chrom_list_file"
    fi
done

echo "染色体文件列表已创建: $chrom_list_file" >> "$log_file"
echo "包含 $(wc -l < "$chrom_list_file") 个染色体文件" >> "$log_file"

# 可选：合并所有染色体的文件（如果需要）
if [ $successful_jobs -eq $chrom_count ] && [ $chrom_count -gt 1 ]; then
    echo "" >> "$log_file"
    echo "所有染色体基因分型成功，可以创建最终合并文件（可选）..." >> "$log_file"
    
    final_merged="$output_dir/genotyped.all_chromosomes.vcf.gz"
    
    echo "合并所有 $chrom_count 条染色体的VCF文件..." >> "$log_file"
    
    $gatk_path --java-options "-Xmx80g -Xms20g" MergeVcfs \
        -R "$reference_genome" \
        --variant "$chrom_list_file" \
        -O "$final_merged" \
        --tmp-dir "$output_dir/tmp/final_merge" >> "$log_file" 2>&1
        
    if [ $? -eq 0 ]; then
        echo "最终合并文件创建成功: $final_merged" >> "$log_file"
        
        # 创建最终索引
        $gatk_path --java-options "-Xmx20g" IndexFeatureFile \
            -I "$final_merged" >> "$log_file" 2>&1
    else
        echo "最终合并文件创建失败" >> "$log_file"
    fi
    
    # 清理临时文件
    rm -rf "$output_dir/tmp/final_merge"
fi

echo "GATK GenotypeGVCFs processing completed at $(date)" >> "$log_file"

# 创建重运行失败的染色体脚本（如果有失败的）
if [ $failed_count -gt 0 ]; then
    retry_script="$output_dir/retry_failed_chromosomes.sh"
    echo "#!/bin/bash" > "$retry_script"
    echo "# 重运行失败的染色体" >> "$retry_script"
    echo "# 生成时间: $(date)" >> "$retry_script"
    echo "" >> "$retry_script"
    
    for chrom in "${chromosomes[@]}"; do
        chrom_file="$output_dir/per_chromosome/genotyped.${chrom}.vcf.gz"
        if [ ! -f "$chrom_file" ] || [ ! -f "${chrom_file}.tbi" ]; then
            echo "echo '=== 重运行染色体 $chrom ==='" >> "$retry_script"
            echo "$GENOTYPE_SCRIPT '$chrom' '$output_dir' '$reference_genome' '$gatk_path' '$input_gvcf_dir'" >> "$retry_script"
            echo "" >> "$retry_script"
        fi
    done
    
    chmod +x "$retry_script"
    echo "已创建重运行失败的染色体脚本: $retry_script" >> "$log_file"
fi

# 在终端显示最终结果
echo ""
echo "=================================================================="
echo "GATK GenotypeGVCFs 按染色体基因分型完成!"
echo "=================================================================="
echo "成功分型: $successful_jobs/${#chromosomes[@]} 条染色体"
echo "详细日志: $log_file"
echo "染色体VCF文件: $output_dir/per_chromosome/"
echo "染色体文件列表: $chrom_list_file"
echo "各染色体日志: $output_dir/logs/"
echo "Parallel任务日志: $output_dir/parallel_jobs.log"

if [ $failed_count -gt 0 ]; then
    echo ""
    echo "注意: 有 $failed_count 条染色体处理失败"
    echo "失败的染色体:"
    grep ": 失败" "$output_dir/genotypegvcfs.log" 2>/dev/null || echo "  无详细记录"
    if [ -f "$retry_script" ]; then
        echo "重运行失败的染色体脚本: $retry_script"
    fi
fi

if [ $successful_jobs -eq $chrom_count ] && [ -f "$output_dir/genotyped.all_chromosomes.vcf.gz" ]; then
    echo ""
    echo "所有染色体合并文件: $output_dir/genotyped.all_chromosomes.vcf.gz"
fi

echo ""
echo "下一步建议:"
echo "1. 使用SelectVariants分离SNP和INDEL"
echo "2. 使用VariantFiltration进行硬过滤"
echo "3. 使用VCFtools或bcftools进行进一步过滤和分析"
echo "=================================================================="
