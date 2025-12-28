#!/bin/bash

# ================= 配置 =================
# gVCF 所在的目录
GVCF_DIR="/home/data/5.gatk_haplotypecaller/gvcf"
# 并发数量 (你服务器32核，考虑到IO压力，设为20比较稳妥且快)
THREADS=20
# 输出文件
OUT_FILE="L_stats.txt"

# ================= 1. 生成列表 =================
echo "正在扫描 gVCF 文件..."
# 后缀 .g.vcf.gz
ls ${GVCF_DIR}/*.g.vcf.gz > full_sample_list.txt

count=$(wc -l < full_sample_list.txt)
echo "共找到 $count 个 gVCF 文件。"

# ================= 2. 定义计算函数 =================
calc_L_single() {
    vcf=$1
    
    # 使用 bcftools query 计算 
    # 逻辑：非变异块长度 = END - POS + 1；变异位点长度 = 1
    L=$(bcftools query -f '%POS\t%END\n' "$vcf" | \
        awk '{
            if ($2 != ".") { len = $2 - $1 + 1 } 
            else { len = 1 }
            sum += len
        } 
        END {printf "%.0f", sum}')
    
    # 输出格式: 文件名 <TAB> 长度
    echo -e "${vcf}\t${L}"
}
export -f calc_L_single

# ================= 3. 并行运行 =================
echo "开始并行计算 L (并发数: $THREADS)..."
echo "结果将写入 $OUT_FILE"

# 使用 GNU Parallel
# --bar: 显示进度条
# --keep-order: 保持输出顺序 (可选)
cat full_sample_list.txt | parallel -j $THREADS --bar calc_L_single {} > "$OUT_FILE"

echo "----------------------------------------"
echo "计算完成！"

# ================= 4. 自动计算平均值 =================
echo "正在计算 135 个样本的平均 L 值..."

# 使用 awk 计算第二列的平均值
AVG_L=$(awk '{sum+=$2} END {printf "%.0f", sum/NR}' "$OUT_FILE")

echo "========================================"
echo "所有样本的 L 统计在: $OUT_FILE"
echo "★ 建议填入 Blueprint 的 L 值 (平均值): $AVG_L"
echo "========================================"
