#!/bin/bash

# 设置目录和文件
input_gvcf_dir="/home/vensin/workspace/snpcalling_wild/5.gatk_haplotypecaller/gvcf"
output_dir="/home/data/6.gatk_combinegvcfs"
reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
gatk_path="/home/vensin/software/gatk-4.6.2.0/gatk"

# 创建输出目录
mkdir -p "$output_dir/per_chromosome"

# 染色体列表
chromosomes=("Chr01" "Chr02" "Chr03" "Chr04" "Chr05" "Chr06" "Chr07" "Chr08" "Chr09" "Chr10"
             "Chr11" "Chr12" "Chr13" "Chr14" "Chr15" "Chr16" "Chr17" "Chr18" "Chr19" "Chr20"
             "Chr21" "Chr22" "Chr23")

# 获取GVCF文件列表
gvcf_files=($(find "$input_gvcf_dir" -name '*.g.vcf.gz'))
gvcf_list_file="$output_dir/gvcf_files.list"
printf "%s\n" "${gvcf_files[@]}" | sort > "$gvcf_list_file"

echo "开始按染色体合并 ${#gvcf_files[@]} 个GVCF文件..."

# 并行处理函数
process_chromosome_parallel() {
    chrom=$1
    echo "合并染色体 $chrom..."
    
    output_file="$output_dir/per_chromosome/combined.${chrom}.g.vcf.gz"
    
    # 合并指定染色体的GVCF
    $gatk_path --java-options "-Xmx20g -Xms10g" CombineGVCFs \
        -R "$reference_genome" \
        --variant "$gvcf_list_file" \
        -L "$chrom" \
        -O "$output_file" \
        --tmp-dir "$output_dir/tmp_${chrom}" \
        --read-index 0
    
    if [ $? -eq 0 ]; then
        # 创建索引
        $gatk_path --java-options "-Xmx4g" IndexFeatureFile \
            -F "$output_file"
        
        echo "染色体 $chrom 合并完成"
        rm -rf "$output_dir/tmp_${chrom}"
    else
        echo "染色体 $chrom 合并失败"
    fi
}

export -f process_chromosome_parallel
export input_gvcf_dir output_dir reference_genome gatk_path gvcf_list_file

# 使用parallel并行处理
if command -v parallel >/dev/null 2>&1; then
    echo "使用GNU Parallel并行处理..."
    
    # 确定并行任务数
    cpu_count=$(nproc)
    parallel_jobs=$(( cpu_count / 2 ))
    if [ $parallel_jobs -lt 1 ]; then
        parallel_jobs=1
    fi
    
    printf "%s\n" "${chromosomes[@]}" | parallel -j $parallel_jobs \
        "process_chromosome_parallel {}"
else
    # 串行处理
    for chrom in "${chromosomes[@]}"; do
        process_chromosome_parallel "$chrom"
    done
fi

# 生成染色体文件列表
chrom_list_file="$output_dir/chromosome_files.list"
> "$chrom_list_file"
for chrom in "${chromosomes[@]}"; do
    chrom_file="$output_dir/per_chromosome/combined.${chrom}.g.vcf.gz"
    if [ -f "$chrom_file" ] && [ -f "${chrom_file}.tbi" ]; then
        echo "$chrom_file" >> "$chrom_list_file"
        file_size=$(ls -lh "$chrom_file" | awk '{print $5}')
        echo "  $chrom: $file_size"
    else
        echo "  $chrom: 文件缺失或索引未创建"
    fi
done

echo "按染色体合并完成!"
echo "染色体合并文件位于: $output_dir/per_chromosome/"
echo "染色体文件列表: $chrom_list_file"
