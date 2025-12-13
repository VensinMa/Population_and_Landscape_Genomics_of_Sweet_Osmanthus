#!/bin/bash

# ======================================================================
# GATK 硬过滤 + 剔除 + 自动清理脚本
# 功能: 
#   1. VariantFiltration (打标签) -> 生成临时 _marked 文件
#   2. SelectVariants (剔除非PASS) -> 生成最终 _final_pass 文件
#   3. 如果第2步成功，自动删除 _marked 中间文件
# ======================================================================

# --- 1. 目录和文件配置 ---
input_base_dir="/home/data/9.gatk_selectvariants"
output_base_dir="/home/data/10.gatk_variant_filtration_pass"

input_snp_vcf="$input_base_dir/SNP/raw_snps.vcf.gz"
input_indel_vcf="$input_base_dir/INDEL/raw_indels.vcf.gz"

snp_out_dir="$output_base_dir/SNP"
indel_out_dir="$output_base_dir/INDEL"

reference_genome="/home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa"
gatk_path="/home/vensin/software/gatk-4.6.2.0/gatk"

# 内存设置 (50G x 2)
JAVA_OPTS="-Xmx50g -Xms10g"

# --- 2. 准备工作 ---
mkdir -p "$output_base_dir" "$snp_out_dir" "$indel_out_dir"
main_log="$output_base_dir/main_process.log"

snp_log="$snp_out_dir/filter_snp_pipeline.log"
indel_log="$indel_out_dir/filter_indel_pipeline.log"

echo "=== GATK 过滤流水线 (含清理) 开始: $(date) ===" | tee "$main_log"

# --- 3. SNP 处理函数 ---
process_snp() {
    echo "=== [SNP] 流程开始 ===" > "$snp_log"
    
    # 临时文件 (带标签)
    marked_snp="$snp_out_dir/snp_marked.vcf.gz"
    # 最终文件 (仅PASS)
    final_snp="$snp_out_dir/snp_filtered.vcf.gz"

    # 步骤 3.1: VariantFiltration
    echo ">> 步骤 1/3: SNP 硬过滤打标签..." >> "$snp_log"
    $gatk_path --java-options "$JAVA_OPTS" VariantFiltration \
        -R "$reference_genome" \
        -V "$input_snp_vcf" \
        -O "$marked_snp" \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \
        --filter-expression "FS > 60.0" --filter-name "FS60" \
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        >> "$snp_log" 2>&1

    if [ $? -ne 0 ]; then echo "[SNP] 打标签失败，停止" >> "$snp_log"; return 1; fi

    # 步骤 3.2: SelectVariants (--exclude-filtered)
    echo ">> 步骤 2/3: 提取 PASS 位点..." >> "$snp_log"
    $gatk_path --java-options "$JAVA_OPTS" SelectVariants \
        -R "$reference_genome" \
        -V "$marked_snp" \
        -O "$final_snp" \
        --exclude-filtered \
        >> "$snp_log" 2>&1
    
    select_exit_code=$?

    # 步骤 3.3: 检查并清理
    if [ $select_exit_code -eq 0 ] && [ -f "$final_snp" ]; then
        echo ">> 步骤 3/3: 最终文件生成成功，删除中间文件..." >> "$snp_log"
        rm -f "$marked_snp" "$marked_snp.tbi"
        echo "   已删除: snp_marked.vcf.gz" >> "$snp_log"
        
        # 统计行数
        count=$(zgrep -v "^#" "$final_snp" | wc -l)
        echo "=== [SNP] 全部完成! 保留位点数: $count ===" >> "$snp_log"
    else
        echo "[SNP] 提取失败，保留中间文件以便排查错误" >> "$snp_log"
        return 1
    fi
}

# --- 4. INDEL 处理函数 ---
process_indel() {
    echo "=== [INDEL] 流程开始 ===" > "$indel_log"
    
    marked_indel="$indel_out_dir/indel_marked.vcf.gz"
    final_indel="$indel_out_dir/indel_filtered.vcf.gz"

    # 步骤 4.1: VariantFiltration
    echo ">> 步骤 1/3: INDEL 硬过滤打标签..." >> "$indel_log"
    $gatk_path --java-options "$JAVA_OPTS" VariantFiltration \
        -R "$reference_genome" \
        -V "$input_indel_vcf" \
        -O "$marked_indel" \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        --filter-expression "FS > 200.0" --filter-name "FS200" \
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        >> "$indel_log" 2>&1

    if [ $? -ne 0 ]; then echo "[INDEL] 打标签失败，停止" >> "$indel_log"; return 1; fi

    # 步骤 4.2: SelectVariants
    echo ">> 步骤 2/3: 提取 PASS 位点..." >> "$indel_log"
    $gatk_path --java-options "$JAVA_OPTS" SelectVariants \
        -R "$reference_genome" \
        -V "$marked_indel" \
        -O "$final_indel" \
        --exclude-filtered \
        >> "$indel_log" 2>&1

    select_exit_code=$?

    # 步骤 4.3: 检查并清理
    if [ $select_exit_code -eq 0 ] && [ -f "$final_indel" ]; then
        echo ">> 步骤 3/3: 最终文件生成成功，删除中间文件..." >> "$indel_log"
        rm -f "$marked_indel" "$marked_indel.tbi"
        echo "   已删除: indel_marked.vcf.gz" >> "$indel_log"
        
        count=$(zgrep -v "^#" "$final_indel" | wc -l)
        echo "=== [INDEL] 全部完成! 保留位点数: $count ===" >> "$indel_log"
    else
        echo "[INDEL] 提取失败，保留中间文件以便排查错误" >> "$indel_log"
        return 1
    fi
}

# --- 5. 并行启动任务 ---
echo ">>> 启动 SNP 处理流程 (后台)..." | tee -a "$main_log"
process_snp &
pid_snp=$!

echo ">>> 启动 INDEL 处理流程 (后台)..." | tee -a "$main_log"
process_indel &
pid_indel=$!

# --- 6. 等待并报告 ---
echo ">>> 正在运行... (完成时将自动清理中间文件)" | tee -a "$main_log"
wait $pid_snp
status_snp=$?
wait $pid_indel
status_indel=$?

echo "" | tee -a "$main_log"
echo "=== 最终结果报告 ===" | tee -a "$main_log"

if [ $status_snp -eq 0 ]; then
    echo "✓ SNP 成功 (中间文件已清理)" | tee -a "$main_log"
    echo "  输出: $snp_out_dir/snp_filtered.vcf.gz" | tee -a "$main_log"
else
    echo "✗ SNP 失败 (中间文件可能保留)" | tee -a "$main_log"
fi

if [ $status_indel -eq 0 ]; then
    echo "✓ INDEL 成功 (中间文件已清理)" | tee -a "$main_log"
    echo "  输出: $indel_out_dir/indel_filtered.vcf.gz" | tee -a "$main_log"
else
    echo "✗ INDEL 失败 (中间文件可能保留)" | tee -a "$main_log"
fi

exit 0
