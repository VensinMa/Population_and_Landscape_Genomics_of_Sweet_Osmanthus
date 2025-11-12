#!/bin/bash

# 设置目录和文件
clean_data_dir="/home/vensin/workspace/snpcalling_wild/2.cleaned_data"
sam_bam_dir="/home/data/3.bwa_sam_bam"
reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"

# 设置 bwa-mem2 路径
BWA_MEM2="/home/vensin/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2"

# 创建输出目录
mkdir -p "$sam_bam_dir/sam" "$sam_bam_dir/bam" "$sam_bam_dir/sorted_bam" "$sam_bam_dir/logs"

# 设置日志文件
log_file="$sam_bam_dir/bwa_processing.log"
echo "BWA-MEM2 processing started at $(date)" > "$log_file"

# 检查样本数量
echo "=== 检查样本数量 ===" >> "$log_file"
sample_count=$(find "$clean_data_dir" -name '*_R1.clean.fastq.gz' | wc -l)
echo "找到 $sample_count 个样本" >> "$log_file"

if [ "$sample_count" -eq 0 ]; then
    echo "错误: 未找到样本文件" | tee -a "$log_file"
    exit 1
fi

# 显示前几个样本
echo "样本示例:" >> "$log_file"
find "$clean_data_dir" -name '*_R1.clean.fastq.gz' | head -5 | xargs -n 1 basename | sed 's/_R1.clean.fastq.gz//' >> "$log_file"

# 构建参考基因组索引（使用 bwa-mem2）
echo "=== 检查参考基因组索引 ===" >> "$log_file"
if [ ! -f "$reference_genome.bwt.2bit.64" ]; then
    echo "构建 BWA-MEM2 索引..." >> "$log_file"
    $BWA_MEM2 index "$reference_genome" >> "$log_file" 2>&1
    if [ $? -ne 0 ]; then
        echo "BWA-MEM2 索引构建失败" | tee -a "$log_file"
        exit 1
    fi
else
    echo "BWA-MEM2 索引已存在" >> "$log_file"
fi

# 其他索引（samtools 和 GATK 的索引仍然需要）
if [ ! -f "$reference_genome.fai" ]; then
    echo "创建 FASTA 索引..." >> "$log_file"
    samtools faidx "$reference_genome" >> "$log_file" 2>&1
    if [ $? -ne 0 ]; then
        echo "FASTA 索引创建失败" | tee -a "$log_file"
        exit 1
    fi
else
    echo "FASTA 索引已存在" >> "$log_file"
fi

dict_file="${reference_genome%.*}.dict"
if [ ! -f "$dict_file" ]; then
    echo "创建序列字典..." >> "$log_file"
    gatk CreateSequenceDictionary -R "$reference_genome" -O "$dict_file" >> "$log_file" 2>&1
    if [ $? -ne 0 ]; then
        echo "序列字典创建失败" | tee -a "$log_file"
        exit 1
    fi
else
    echo "序列字典已存在" >> "$log_file"
fi

# 定义处理函数（使用 bwa-mem2）
process_sample() {
    local fq1=$1
    local base_name=$(basename "$fq1" "_R1.clean.fastq.gz")
    local fq2="${fq1/_R1.clean.fastq.gz/_R2.clean.fastq.gz}"
    
    # 样本特定的日志
    local sample_log="$sam_bam_dir/logs/${base_name}.log"
    
    echo "开始处理样本: $base_name" >> "$sample_log"
    echo "  R1: $(basename "$fq1")" >> "$sample_log"
    echo "  R2: $(basename "$fq2")" >> "$sample_log"
    
    # 检查R2文件是否存在
    if [ ! -f "$fq2" ]; then
        echo "错误: 找不到对应的R2文件: $fq2" >> "$sample_log"
        return 1
    fi
    
    # 执行 BWA-MEM2 比对
    echo "  BWA-MEM2 比对..." >> "$sample_log"
    if ! $BWA_MEM2 mem -t 4 -M -R "@RG\tID:$base_name\tSM:$base_name\tLB:$base_name\tPL:ILLUMINA" "$reference_genome" "$fq1" "$fq2" > "$sam_bam_dir/sam/${base_name}.sam" 2>> "$sample_log"; then
        echo "  BWA-MEM2 失败: $base_name" >> "$sample_log"
        return 1
    fi
    
    # 转换SAM为BAM
    echo "  SAM转BAM..." >> "$sample_log"
    if ! samtools view -@ 4 -bS "$sam_bam_dir/sam/${base_name}.sam" -o "$sam_bam_dir/bam/${base_name}.bam" 2>> "$sample_log"; then
        echo "  SAM转BAM失败: $base_name" >> "$sample_log"
        return 1
    fi
    
    # 删除SAM文件以节省空间
    rm "$sam_bam_dir/sam/${base_name}.sam"
    
    # 排序BAM文件
    echo "  BAM排序..." >> "$sample_log"
    if ! samtools sort -@ 4 -m 15G "$sam_bam_dir/bam/${base_name}.bam" -o "$sam_bam_dir/sorted_bam/${base_name}.sorted.bam" 2>> "$sample_log"; then
        echo "  BAM排序失败: $base_name" >> "$sample_log"
        return 1
    fi
    
    # 删除未排序的BAM文件
    rm "$sam_bam_dir/bam/${base_name}.bam"
    
    # 为排序后的BAM文件创建索引
    echo "  BAM索引..." >> "$sample_log"
    if ! samtools index "$sam_bam_dir/sorted_bam/${base_name}.sorted.bam" 2>> "$sample_log"; then
        echo "  BAM索引失败: $base_name" >> "$sample_log"
        return 1
    fi
    
    echo "  完成: $base_name" >> "$sample_log"
    echo "$base_name: 成功" >> "$log_file"
    return 0
}

export -f process_sample
export clean_data_dir
export sam_bam_dir
export reference_genome
export BWA_MEM2

# 处理所有样本
echo "=== 开始处理样本 ===" >> "$log_file"

# 使用parallel处理（根据系统资源调整并行数）
find "$clean_data_dir" -name '*_R1.clean.fastq.gz' | sort | \
    parallel --joblog "$sam_bam_dir/parallel_jobs.log" -j 8 \
    "process_sample {}"

# 检查处理结果
successful_jobs=$(grep -c ": 成功" "$log_file" 2>/dev/null || echo 0)
echo "处理完成: $successful_jobs/$sample_count 个样本成功处理" >> "$log_file"

# 显示处理摘要
echo "=== 处理摘要 ===" >> "$log_file"
echo "成功: $successful_jobs" >> "$log_file"
echo "失败: $((sample_count - successful_jobs))" >> "$log_file"

echo "BWA-MEM2 processing completed at $(date)" >> "$log_file"

# 在终端显示最终结果
echo "处理完成!"
echo "成功处理: $successful_jobs/$sample_count 个样本"
echo "详细日志: $log_file"
echo "各样本日志: $sam_bam_dir/logs/"
