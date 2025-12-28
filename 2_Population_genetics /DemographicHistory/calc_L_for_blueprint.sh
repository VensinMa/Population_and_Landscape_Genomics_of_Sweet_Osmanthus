#!/bin/bash

# ================= 配置 =================
# gVCF 所在的目录
GVCF_DIR="/home/data/5.gatk_haplotypecaller/gvcf"

# 并发数量
THREADS=20

# 输出结果文件
OUT_FILE="L_stats_selected.txt"

# 手动指定的样本ID列表文件 (单列，一行一个样本ID，无表头行名等)
# 如果文件不存在，脚本会报错提示
MANUAL_ID_FILE="sample_id_list.txt"

# gVCF 文件的后缀 (用于拼接路径)
SUFFIX=".g.vcf.gz"

# 临时生成的待计算文件列表
TARGET_LIST="target_files.txt"

# ================= 1. 生成待计算列表 =================
echo "正在生成待计算文件列表..."
> "$TARGET_LIST" # 清空旧列表

# --- 模式 A：使用手动指定样本 ID (当前启用) ---
if [ -f "$MANUAL_ID_FILE" ]; then
    echo "检测到样本列表文件: $MANUAL_ID_FILE"
    while read -r sample_id; do
        # 跳过空行
        [[ -z "$sample_id" ]] && continue
        
        # 拼接完整路径: 目录 + 样本ID + 后缀
        full_path="${GVCF_DIR}/${sample_id}${SUFFIX}"
        
        if [ -f "$full_path" ]; then
            echo "$full_path" >> "$TARGET_LIST"
        else
            echo "警告: 找不到文件 $full_path ，已跳过。"
        fi
    done < "$MANUAL_ID_FILE"
else
    echo "错误: 未找到样本列表文件 $MANUAL_ID_FILE"
    echo "请创建该文件，或取消下方'模式 B'的注释以计算目录下所有文件。"
    exit 1
fi

# --- 模式 B：计算目录下所有 gVCF (已注释，如需使用请取消注释并注释掉模式 A) ---
# echo "正在扫描目录下所有 gVCF 文件..."
# ls ${GVCF_DIR}/*${SUFFIX} > "$TARGET_LIST"


# 检查列表是否为空
count=$(wc -l < "$TARGET_LIST")
if [ "$count" -eq 0 ]; then
    echo "错误: 待计算列表为空！请检查样本ID是否正确或文件是否存在。"
    exit 1
fi
echo "共找到 $count 个有效 gVCF 文件准备计算。"


# ================= 2. 定义计算函数  =================
calc_L_single() {
    vcf=$1
    
    # 使用 bcftools query 计算 
    # 逻辑：非变异块长度 = END - POS + 1；变异位点长度 = 1
    # 注意：这里去掉了 --threads 以兼容不同版本
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
cat "$TARGET_LIST" | parallel -j $THREADS --bar calc_L_single {} > "$OUT_FILE"

echo "----------------------------------------"
echo "计算完成！"

# ================= 4. 自动计算平均值 =================
echo "正在计算所选样本的平均 L 值..."

# 使用 awk 计算第二列的平均值
AVG_L=$(awk '{sum+=$2} END {printf "%.0f", sum/NR}' "$OUT_FILE")

echo "========================================"
echo "统计样本数: $count"
echo "详细结果文件: $OUT_FILE"
echo "★ 建议填入 Blueprint 的 L 值 (平均值): $AVG_L"
echo "========================================"
