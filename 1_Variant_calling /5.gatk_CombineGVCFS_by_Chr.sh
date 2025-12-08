#!/bin/bash

# 设置目录和文件
input_gvcf_dir="/home/vensin/workspace/snpcalling_wild/5.gatk_haplotypecaller/gvcf"
output_dir="/home/data/6.gatk_combinegvcf_per_chr"
reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
gatk_path="/home/vensin/software/gatk-4.6.2.0/gatk"

# 默认并行任务数（可根据需要修改）
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
mkdir -p "$output_dir/per_chromosome" "$output_dir/logs" "$output_dir/tmp"

# 染色体列表
chromosomes=("Chr01" "Chr02" "Chr03" "Chr04" "Chr05" "Chr06" "Chr07" "Chr08" "Chr09" "Chr10"
             "Chr11" "Chr12" "Chr13" "Chr14" "Chr15" "Chr16" "Chr17" "Chr18" "Chr19" "Chr20"
             "Chr21" "Chr22" "Chr23")

# 获取GVCF文件列表
gvcf_files=($(find "$input_gvcf_dir" -name '*.g.vcf.gz'))
gvcf_count=${#gvcf_files[@]}

if [ "$gvcf_count" -eq 0 ]; then
    echo "错误: 未找到GVCF文件"
    exit 1
fi

echo "找到 $gvcf_count 个GVCF文件"

gvcf_list_file="$output_dir/gvcf_files.list"
printf "%s\n" "${gvcf_files[@]}" | sort > "$gvcf_list_file"

echo "开始按染色体合并 $gvcf_count 个GVCF文件..."

# 并行处理函数
process_chromosome() {
    chrom=$1
    task_num=$2
    
    # 为每个任务创建独立日志文件
    chrom_log="$output_dir/logs/combine_${chrom}.log"
    echo "=== 处理染色体 $chrom (任务 $task_num) ===" > "$chrom_log"
    echo "开始时间: $(date)" >> "$chrom_log"
    
    echo "[任务 $task_num] 合并染色体 $chrom..."
    
    output_file="$output_dir/per_chromosome/combined.${chrom}.g.vcf.gz"
    
    # 创建染色体特定的临时目录
    chrom_tmp_dir="$output_dir/tmp/${chrom}"
    mkdir -p "$chrom_tmp_dir"
    
    echo "合并指定染色体的GVCF..." >> "$chrom_log"
    
    # 合并指定染色体的GVCF
    $gatk_path --java-options "-Xmx40g -Xms20g" CombineGVCFs \
        -R "$reference_genome" \
        --variant "$gvcf_list_file" \
        -L "$chrom" \
        -O "$output_file" \
        --tmp-dir "$chrom_tmp_dir" \
        --read-index 0 >> "$chrom_log" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "[任务 $task_num] 染色体 $chrom 合并完成，开始创建索引..." >> "$chrom_log"
        
        # 创建索引
        $gatk_path --java-options "-Xmx20g -Xms10g" IndexFeatureFile \
            -F "$output_file" >> "$chrom_log" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "[任务 $task_num] 染色体 $chrom 索引创建成功" >> "$chrom_log"
            echo "染色体 $chrom 合并完成"
            
            # 清理临时目录
            rm -rf "$chrom_tmp_dir"
        else
            echo "[任务 $task_num] 染色体 $chrom 索引创建失败" >> "$chrom_log"
            echo "染色体 $chrom 索引创建失败"
            return 1
        fi
    else
        echo "[任务 $task_num] 染色体 $chrom 合并失败" >> "$chrom_log"
        echo "染色体 $chrom 合并失败"
        return 1
    fi
    
    echo "结束时间: $(date)" >> "$chrom_log"
    return 0
}

export -f process_chromosome
export input_gvcf_dir output_dir reference_genome gatk_path gvcf_list_file

# 使用parallel并行处理
if command -v parallel >/dev/null 2>&1; then
    echo "使用 GNU Parallel 并行处理染色体..."
    
    # 确定并行任务数（不超过染色体数）
    max_jobs=${#chromosomes[@]}
    if [ "$PARALLEL_JOBS" -gt "$max_jobs" ]; then
        echo "警告: 并行任务数($PARALLEL_JOBS)超过染色体数($max_jobs)，调整为 $max_jobs"
        PARALLEL_JOBS="$max_jobs"
    fi
    
    echo "并行处理 ${#chromosomes[@]} 条染色体，并行任务数: $PARALLEL_JOBS"
    
    # 为每个染色体分配任务编号
    task_index=0
    chrom_tasks=()
    for chrom in "${chromosomes[@]}"; do
        ((task_index++))
        chrom_tasks+=("$chrom $task_index")
    done
    
    # 并行处理所有染色体
    printf "%s\n" "${chrom_tasks[@]}" | parallel -j $PARALLEL_JOBS \
        --colsep ' ' \
        --joblog "$output_dir/chromosome_jobs.log" \
        --progress \
        --eta \
        "process_chromosome {1} {2}"
    
    parallel_exit_code=$?
    
    if [ $parallel_exit_code -ne 0 ]; then
        echo "部分染色体处理失败，请查看日志文件"
    fi
else
    # 串行处理（不使用parallel）
    echo "GNU Parallel 未安装，使用串行处理"
    echo "注意: 串行处理较慢，建议安装 GNU Parallel"
    
    task_index=0
    for chrom in "${chromosomes[@]}"; do
        ((task_index++))
        process_chromosome "$chrom" "$task_index"
    done
fi

# 统计成功和失败的染色体
success_count=0
fail_count=0
failed_chromosomes=()

for chrom in "${chromosomes[@]}"; do
    chrom_file="$output_dir/per_chromosome/combined.${chrom}.g.vcf.gz"
    if [ -f "$chrom_file" ] && [ -f "${chrom_file}.tbi" ]; then
        ((success_count++))
    else
        ((fail_count++))
        failed_chromosomes+=("$chrom")
    fi
done

# 生成染色体文件列表
chrom_list_file="$output_dir/chromosome_files.list"
> "$chrom_list_file"

echo ""
echo "=== 合并结果统计 ==="
echo "成功合并的染色体: $success_count/${#chromosomes[@]}"
echo "失败的染色体: $fail_count"

if [ $fail_count -gt 0 ]; then
    echo "失败的染色体列表:"
    for chrom in "${failed_chromosomes[@]}"; do
        echo "  - $chrom"
        echo "    日志文件: $output_dir/logs/combine_${chrom}.log"
    done
fi

echo ""
echo "生成的染色体文件:"

for chrom in "${chromosomes[@]}"; do
    chrom_file="$output_dir/per_chromosome/combined.${chrom}.g.vcf.gz"
    if [ -f "$chrom_file" ] && [ -f "${chrom_file}.tbi" ]; then
        file_size=$(ls -lh "$chrom_file" 2>/dev/null | awk '{print $5}' || echo "未知")
        echo "  $chrom: $file_size"
        echo "$chrom_file" >> "$chrom_list_file"
    else
        echo "  $chrom: 文件缺失或索引未创建"
    fi
done

echo ""
echo "按染色体合并完成!"
echo "染色体合并文件位于: $output_dir/per_chromosome/"
echo "染色体文件列表: $chrom_list_file"
echo "日志文件目录: $output_dir/logs/"
echo "任务日志: $output_dir/chromosome_jobs.log"

# 如果所有染色体都成功合并，创建合并文件列表
if [ $success_count -eq ${#chromosomes[@]} ]; then
    echo "所有染色体合并成功！"
    echo "可以在后续分析中使用染色体文件列表: $chrom_list_file"
elif [ $fail_count -gt 0 ]; then
    echo "有 $fail_count 条染色体合并失败，请检查日志文件重新运行失败的染色体"
    
    # 创建失败染色体的重运行脚本
    if [ ${#failed_chromosomes[@]} -gt 0 ]; then
        retry_script="$output_dir/retry_failed_chromosomes.sh"
        echo "#!/bin/bash" > "$retry_script"
        echo "# 重运行失败的染色体" >> "$retry_script"
        echo "" >> "$retry_script"
        
        for chrom in "${failed_chromosomes[@]}"; do
            echo "echo '重运行染色体 $chrom...'" >> "$retry_script"
            echo "process_chromosome '$chrom' 'retry'" >> "$retry_script"
        done
        
        chmod +x "$retry_script"
        echo ""
        echo "已创建重运行脚本: $retry_script"
        echo "用法: $retry_script"
    fi
fi
