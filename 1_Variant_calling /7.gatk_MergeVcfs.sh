#!/bin/bash

# ======================================================================
# GATK MergeVcfs 合并染色体VCF文件脚本
# 作者：基于上一步GenotypeGVCFs脚本修改
# ======================================================================

# 设置目录和文件
input_vcf_dir="/home/data/7.gatk_genotypegvcfs_per_chr/per_chromosome"  # 上一步按染色体基因分型的VCF文件
output_dir="/home/data/8.gatk_mergevcfs"
reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
gatk_path="/home/vensin/software/gatk-4.6.2.0/gatk"

# 默认并行任务数（仅用于可选的分步处理）
PARALLEL_JOBS=1  # MergeVcfs通常是单任务，但保留参数用于兼容性

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -j|--jobs)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        -f|--force)
            FORCE_RUN=true
            shift
            ;;
        -h|--help)
            echo "用法: $0 [-j 并行任务数] [-f]"
            echo "   -j, --jobs     指定并行任务数（默认: 1，仅用于兼容性）"
            echo "   -f, --force    强制重新运行，即使输出文件已存在"
            echo ""
            echo "输入文件目录: $input_vcf_dir"
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

echo "开始合并染色体VCF文件..."

# 创建输出目录
mkdir -p "$output_dir" "$output_dir/logs"

# 设置日志文件
log_file="$output_dir/mergevcfs.log"
echo "GATK MergeVcfs 合并开始时间: $(date)" > "$log_file"

# 检查输入VCF文件（上一步按染色体基因分型的VCF文件）
echo "=== 检查输入VCF文件 ===" >> "$log_file"
vcf_count=$(find "$input_vcf_dir" -name 'genotyped.*.vcf.gz' | wc -l)

if [ "$vcf_count" -eq 0 ]; then
    echo "错误: 未找到按染色体基因分型的VCF文件" | tee -a "$log_file"
    echo "请检查输入目录: $input_vcf_dir" | tee -a "$log_file"
    echo "预期文件格式: genotyped.Chr01.vcf.gz, genotyped.Chr02.vcf.gz 等" | tee -a "$log_file"
    exit 1
fi

echo "找到 $vcf_count 个按染色体基因分型的VCF文件" >> "$log_file"

# 显示前几个染色体文件
echo "染色体VCF文件示例:" >> "$log_file"
find "$input_vcf_dir" -name 'genotyped.*.vcf.gz' | head -10 | xargs -n 1 basename >> "$log_file"

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
for vcf_file in "$input_vcf_dir"/genotyped.*.vcf.gz; do
    if [ -f "$vcf_file" ]; then
        # 提取染色体名，例如从 genotyped.Chr01.vcf.gz 提取 Chr01
        chrom=$(basename "$vcf_file" .vcf.gz | sed 's/genotyped\.//')
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
IFS=$'\n' chromosomes=($(sort <<<"${chromosomes[@]}"))
unset IFS

echo "找到 ${#chromosomes[@]} 条染色体需要合并:" >> "$log_file"
printf "%s\n" "${chromosomes[@]}" >> "$log_file"

# 创建染色体文件列表（用于合并）
chrom_list_file="$output_dir/chromosome_files.list"
> "$chrom_list_file"
for chrom in "${chromosomes[@]}"; do
    chrom_file="$input_vcf_dir/genotyped.${chrom}.vcf.gz"
    if [ -f "$chrom_file" ] && [ -f "${chrom_file}.tbi" ]; then
        echo "$chrom_file" >> "$chrom_list_file"
        echo "  已添加: $chrom_file" >> "$log_file"
    else
        echo "警告: 染色体 $chrom 的VCF文件或索引文件不存在" >> "$log_file"
        echo "  文件: $chrom_file" >> "$log_file"
        echo "  索引: ${chrom_file}.tbi" >> "$log_file"
    fi
done

# 统计可合并的文件数量
mergeable_count=$(wc -l < "$chrom_list_file")
echo "可合并的文件数量: $mergeable_count/${#chromosomes[@]}" >> "$log_file"

if [ "$mergeable_count" -eq 0 ]; then
    echo "错误: 没有找到可合并的VCF文件（所有文件都缺少或索引不完整）" | tee -a "$log_file"
    exit 1
fi

if [ "$mergeable_count" -ne "${#chromosomes[@]}" ]; then
    echo "警告: 有 $((chrom_count - mergeable_count)) 条染色体的文件无法合并" >> "$log_file"
    echo "继续合并可用的 $mergeable_count 条染色体..." >> "$log_file"
fi

# 输出文件
final_merged="$output_dir/merged.all_chromosomes.vcf.gz"

# 检查输出文件是否已存在
if [ -f "$final_merged" ] && [ -f "${final_merged}.tbi" ] && [ "${FORCE_RUN}" != "true" ]; then
    echo "输出文件已存在: $final_merged" | tee -a "$log_file"
    echo "使用 -f 选项强制重新运行" | tee -a "$log_file"
    
    # 显示现有文件信息
    file_size=$(ls -lh "$final_merged" | awk '{print $5}')
    echo "现有合并文件大小: $file_size" >> "$log_file"
    echo "跳过合并步骤" >> "$log_file"
    
    echo "GATK MergeVcfs 合并完成时间: $(date)" >> "$log_file"
    
    # 在终端显示最终结果
    echo ""
    echo "=================================================================="
    echo "GATK MergeVcfs 合并染色体VCF文件完成!"
    echo "=================================================================="
    echo "合并文件已存在: $final_merged"
    echo "文件大小: $file_size"
    echo "详细日志: $log_file"
    echo ""
    echo "下一步建议:"
    echo "1. 使用SelectVariants分离SNP和INDEL"
    echo "2. 使用VariantFiltration进行硬过滤"
    echo "3. 使用VCFtools或bcftools进行进一步过滤和分析"
    echo "=================================================================="
    exit 0
fi

# 执行 GATK MergeVcfs 合并所有染色体
echo "=== 开始GATK MergeVcfs合并染色体VCF文件 ===" >> "$log_file"
echo "合并 $mergeable_count 条染色体的VCF文件..." >> "$log_file"

# 创建临时目录
tmp_dir="$output_dir/tmp"
mkdir -p "$tmp_dir"

# 检查是否有足够的磁盘空间
echo "检查磁盘空间..." >> "$log_file"
available_space=$(df "$output_dir" | awk 'NR==2 {print $4}')
if [ "$available_space" -lt 10000000 ]; then  # 小于10GB
    echo "警告: 可用磁盘空间不足 (少于10GB)" >> "$log_file"
    echo "可用空间: $available_space KB" >> "$log_file"
else
    echo "可用磁盘空间: $available_space KB" >> "$log_file"
fi

# 执行合并
echo "运行MergeVcfs..." >> "$log_file"
$gatk_path --java-options "-Xmx80g -Xms20g" MergeVcfs \
    -R "$reference_genome" \
    --variant "$chrom_list_file" \
    -O "$final_merged" \
    --tmp-dir "$tmp_dir" >> "$log_file" 2>&1

merge_exit_code=$?

if [ $merge_exit_code -eq 0 ]; then
    echo "MergeVcfs 合并成功!" >> "$log_file"
    
    # 验证输出的VCF文件
    if [ -f "$final_merged" ]; then
        # 获取合并后的文件大小
        file_size=$(ls -lh "$final_merged" | awk '{print $5}')
        echo "合并文件大小: $file_size" >> "$log_file"
        
        # 创建VCF索引
        echo "创建合并VCF索引..." >> "$log_file"
        $gatk_path --java-options "-Xmx20g" IndexFeatureFile \
            -I "$final_merged" >> "$log_file" 2>&1
            
        if [ $? -eq 0 ]; then
            echo "索引创建成功" >> "$log_file"
            
            # 验证合并文件
            echo "验证合并文件..." >> "$log_file"
            $gatk_path --java-options "-Xmx20g" ValidateVariants \
                -R "$reference_genome" \
                -V "$final_merged" >> "$log_file" 2>&1 || true
        else
            echo "索引创建失败" >> "$log_file"
            merge_exit_code=1
        fi
    else
        echo "合并VCF文件未生成" >> "$log_file"
        merge_exit_code=1
    fi
else
    echo "MergeVcfs 合并失败 (退出码: $merge_exit_code)" >> "$log_file"
fi

# 清理临时文件
rm -rf "$tmp_dir"
echo "GATK MergeVcfs processing completed at $(date)" >> "$log_file"

# 创建重运行脚本（如果需要）
if [ $merge_exit_code -ne 0 ]; then
    retry_script="$output_dir/retry_mergevcfs.sh"
    echo "#!/bin/bash" > "$retry_script"
    echo "# 重运行合并VCF文件" >> "$retry_script"
    echo "# 生成时间: $(date)" >> "$retry_script"
    echo "" >> "$retry_script"
    echo "echo '=== 重新运行 MergeVcfs ==='" >> "$retry_script"
    echo "$0 -f" >> "$retry_script"
    chmod +x "$retry_script"
    echo "已创建重运行脚本: $retry_script" >> "$log_file"
fi

# 在终端显示最终结果
echo ""
echo "=================================================================="
echo "GATK MergeVcfs 合并染色体VCF文件完成!"
echo "=================================================================="

if [ $merge_exit_code -eq 0 ]; then
    echo "✓ 成功合并 $mergeable_count 条染色体的VCF文件"
    echo "合并文件: $final_merged"
    
    if [ -f "$final_merged" ]; then
        file_size=$(ls -lh "$final_merged" | awk '{print $5}')
        echo "文件大小: $file_size"
        
        if [ -n "$variant_count" ] && [ "$variant_count" != "未知" ]; then
            echo "总变异位点数: $variant_count"
        fi
    fi
else
    echo "✗ 合并失败"
    echo "退出码: $merge_exit_code"
    if [ -f "$retry_script" ]; then
        echo "重运行脚本: $retry_script"
    fi
fi

echo "详细日志: $log_file"
echo "输入文件列表: $chrom_list_file"
echo "输入染色体: ${#chromosomes[@]} 条"
echo "成功合并: $mergeable_count 条"

echo ""
echo "下一步建议:"
echo "1. 使用SelectVariants分离SNP和INDEL"
echo "2. 使用VariantFiltration进行硬过滤"
echo "3. 使用VCFtools或bcftools进行进一步过滤和分析"
echo "=================================================================="

exit $merge_exit_code
