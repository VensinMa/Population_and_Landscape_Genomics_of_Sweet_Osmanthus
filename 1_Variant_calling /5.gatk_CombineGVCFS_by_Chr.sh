#!/bin/bash

# 设置目录和文件
input_gvcf_dir="/home/vensin/workspace/snpcalling_wild/5.gatk_haplotypecaller/gvcf"
output_dir="/home/data/6.gatk_combinegvcf_per_chr"
reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
gatk_path="/home/vensin/software/gatk-4.6.2.0/gatk"

# 默认并行任务数
PARALLEL_JOBS=4

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -j|--jobs)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        -h|--help)
            echo "用法: $0 [-j 并行任务数]"
            echo "   -j, --jobs     指定并行任务数（默认: 4）"
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
mkdir -p "$output_dir" "$output_dir/per_chromosome" "$output_dir/logs"

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

# 创建GVCF文件列表
gvcf_list_file="$output_dir/gvcf_files.list"
printf "%s\n" "${gvcf_files[@]}" | sort > "$gvcf_list_file"

# 染色体列表
chromosomes=("Chr01" "Chr02" "Chr03" "Chr04" "Chr05" "Chr06" "Chr07" "Chr08" "Chr09" "Chr10"
             "Chr11" "Chr12" "Chr13" "Chr14" "Chr15" "Chr16" "Chr17" "Chr18" "Chr19" "Chr20"
             "Chr21" "Chr22" "Chr23")

echo "开始按染色体合并 $gvcf_count 个GVCF文件..." >> "$log_file"
echo "共 ${#chromosomes[@]} 条染色体需要处理" >> "$log_file"

# 染色体处理函数
process_chromosome() {
    local chrom="$1"
    local output_file="$output_dir/per_chromosome/combined.${chrom}.g.vcf.gz"
    local chrom_log="$output_dir/logs/combine_${chrom}.log"
    local chrom_tmp_dir="$output_dir/tmp_${chrom}"
    
    echo "处理染色体 $chrom 开始时间: $(date)" > "$chrom_log"
    echo "合并染色体 $chrom..." >> "$chrom_log"
    
    mkdir -p "$chrom_tmp_dir"
    
    # 执行CombineGVCFs（按染色体）
    $gatk_path --java-options "-Xmx40g -Xms20g" CombineGVCFs \
        -R "$reference_genome" \
        --variant "$gvcf_list_file" \
        -L "$chrom" \
        -O "$output_file" \
        --tmp-dir "$chrom_tmp_dir" \
        --read-index 0 >> "$chrom_log" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "染色体 $chrom CombineGVCFs 成功完成" >> "$chrom_log"
        
        # 创建GVCF索引
        echo "创建染色体 $chrom 的GVCF索引..." >> "$chrom_log"
        $gatk_path --java-options "-Xmx20g -Xms10g" IndexFeatureFile \
            -F "$output_file" >> "$chrom_log" 2>&1
            
        if [ $? -eq 0 ]; then
            echo "染色体 $chrom 索引创建完成" >> "$chrom_log"
            
            # 基本验证
            echo "验证染色体 $chrom 的合并文件..." >> "$chrom_log"
            $gatk_path --java-options "-Xmx10g" ValidateVariants \
                -V "$output_file" \
                -R "$reference_genome" \
                -L "$chrom" >> "$chrom_log" 2>&1
                
            if [ $? -eq 0 ]; then
                echo "染色体 $chrom GVCF文件验证通过" >> "$chrom_log"
            else
                echo "警告: 染色体 $chrom GVCF文件验证失败，但文件可能仍可用" >> "$chrom_log"
            fi
            
            # 收集统计信息
            echo "收集染色体 $chrom 的基本统计信息..." >> "$chrom_log"
            variant_count=$($gatk_path --java-options "-Xmx10g" CountVariants \
                -V "$output_file" \
                -L "$chrom" 2>/dev/null | head -1)
            
            echo "染色体 $chrom 统计:" >> "$chrom_log"
            echo "样本数量: $gvcf_count" >> "$chrom_log"
            echo "变异位点数: ${variant_count:-未知}" >> "$chrom_log"
            echo "输出文件: $output_file" >> "$chrom_log"
        else
            echo "染色体 $chrom 索引创建失败" >> "$chrom_log"
            echo "染色体 $chrom 处理失败" >> "$log_file"
            rm -rf "$chrom_tmp_dir"
            return 1
        fi
    else
        echo "染色体 $chrom CombineGVCFs 失败" >> "$chrom_log"
        echo "染色体 $chrom 处理失败" >> "$log_file"
        rm -rf "$chrom_tmp_dir"
        return 1
    fi
    
    # 清理临时目录
    rm -rf "$chrom_tmp_dir"
    
    echo "染色体 $chrom 处理完成时间: $(date)" >> "$chrom_log"
    echo "染色体 $chrom 处理成功" >> "$log_file"
    return 0
}

export -f process_chromosome
export input_gvcf_dir output_dir reference_genome gatk_path gvcf_list_file gvcf_count log_file

# 使用GNU Parallel并行处理
if command -v parallel >/dev/null 2>&1; then
    echo "使用 GNU Parallel 并行处理染色体..." >> "$log_file"
    
    # 确定并行任务数（不超过染色体数）
    max_jobs=${#chromosomes[@]}
    if [ "$PARALLEL_JOBS" -gt "$max_jobs" ]; then
        echo "并行任务数($PARALLEL_JOBS)超过染色体数($max_jobs)，调整为 $max_jobs" >> "$log_file"
        PARALLEL_JOBS="$max_jobs"
    fi
    
    echo "并行处理 ${#chromosomes[@]} 条染色体，并行任务数: $PARALLEL_JOBS" >> "$log_file"
    
    # 并行处理所有染色体
    printf "%s\n" "${chromosomes[@]}" | parallel -j $PARALLEL_JOBS \
        --joblog "$output_dir/chromosome_jobs.log" \
        --progress \
        "process_chromosome {}"
    
    parallel_exit_code=$?
    
    if [ $parallel_exit_code -ne 0 ]; then
        echo "部分染色体处理失败" >> "$log_file"
    fi
else
    # 串行处理
    echo "GNU Parallel 未安装，使用串行处理..." >> "$log_file"
    for chrom in "${chromosomes[@]}"; do
        echo "处理染色体 $chrom..." >> "$log_file"
        process_chromosome "$chrom"
    done
fi

# 统计成功和失败的染色体
success_count=0
fail_count=0
for chrom in "${chromosomes[@]}"; do
    chrom_file="$output_dir/per_chromosome/combined.${chrom}.g.vcf.gz"
    if [ -f "$chrom_file" ] && [ -f "${chrom_file}.tbi" ]; then
        ((success_count++))
    else
        ((fail_count++))
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

# 最终统计
echo "按染色体合并结果统计:" >> "$log_file"
echo "成功合并的染色体数: $success_count/${#chromosomes[@]}" >> "$log_file"
echo "失败的染色体数: $fail_count" >> "$log_file"

echo "按染色体合并的GVCF文件位于: $output_dir/per_chromosome/" >> "$log_file"
echo "染色体文件列表: $chrom_list_file" >> "$log_file"

# 如果有失败的染色体，创建重运行脚本
if [ $fail_count -gt 0 ]; then
    retry_script="$output_dir/retry_failed_chromosomes.sh"
    echo "#!/bin/bash" > "$retry_script"
    echo "# 重运行失败的染色体" >> "$retry_script"
    echo "" >> "$retry_script"
    
    for chrom in "${chromosomes[@]}"; do
        chrom_file="$output_dir/per_chromosome/combined.${chrom}.g.vcf.gz"
        if [ ! -f "$chrom_file" ] || [ ! -f "${chrom_file}.tbi" ]; then
            echo "echo '重运行染色体 $chrom...'" >> "$retry_script"
            echo "process_chromosome '$chrom'" >> "$retry_script"
        fi
    done
    
    chmod +x "$retry_script"
    echo "已创建重运行失败的染色体脚本: $retry_script" >> "$log_file"
fi

# 可选：合并所有染色体的文件
if [ $success_count -eq ${#chromosomes[@]} ]; then
    echo "所有染色体合并成功，开始合并所有染色体文件..." >> "$log_file"
    
    final_combined="$output_dir/combined.all_chromosomes.g.vcf.gz"
    
    $gatk_path --java-options "-Xmx80g -Xms40g" MergeVcfs \
        -R "$reference_genome" \
        --variant "$chrom_list_file" \
        -O "$final_combined" \
        --tmp-dir "$output_dir/tmp_final_merge" >> "$log_file" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "最终合并文件创建成功: $final_combined" >> "$log_file"
        
        # 创建最终索引
        $gatk_path --java-options "-Xmx20g" IndexFeatureFile \
            -F "$final_combined" >> "$log_file" 2>&1
    else
        echo "最终合并文件创建失败" >> "$log_file"
    fi
    
    # 清理临时目录
    rm -rf "$output_dir/tmp_final_merge"
fi

echo "GATK CombineGVCFs 按染色体合并完成时间: $(date)" >> "$log_file"

# 显示结果
echo ""
echo "GATK CombineGVCFs 按染色体合并完成!"
echo "成功合并 $success_count/${#chromosomes[@]} 条染色体"
echo "染色体合并文件位于: $output_dir/per_chromosome/"
echo "染色体文件列表: $chrom_list_file"
echo "日志文件: $log_file"

if [ -f "$retry_script" ]; then
    echo "重运行失败的染色体脚本: $retry_script"
fi

if [ $fail_count -gt 0 ]; then
    echo "注意: 有 $fail_count 条染色体合并失败，请检查日志文件"
    echo "失败的染色体日志文件: $output_dir/logs/combine_*.log"
fi
