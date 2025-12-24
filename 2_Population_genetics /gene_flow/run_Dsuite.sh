#!/bin/bash
set -e  # 遇到错误立即停止

# ==========================================
# 1. 变量设置与目录准备
# ==========================================
WORK_DIR="/home/vensin/workspace/snpcalling_wild/14.gene_flow/Dsuite"
INPUT_VCF="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/208_samples_snp_filtered.LD.pruned.recode.vcf.gz"
POP_FILE_SOURCE="/home/vensin/workspace/snpcalling_wild/14.gene_flow/Dsuite/208samples.pop"

# 软件路径
DSUITE_BIN="/home/vensin/software/Dsuite/Build/Dsuite"
DTOOLS_PY="/home/vensin/software/Dsuite/utils/dtools.py"

echo "=== 1. 初始化工作目录 ==="
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# 检查源文件是否存在
if [ ! -f "$POP_FILE_SOURCE" ]; then
    echo "错误：找不到群体映射文件 $POP_FILE_SOURCE"
    exit 1
fi

# ==========================================
# 2. 生成配置文件 (SETS.txt)
# ==========================================
echo "=== 2. 生成 SETS 配置文件 ==="

# 生成群体水平配置 (SETS_pop.txt): 第1列样本, 第2列群体
# 将 OUTGROUP 替换为 Outgroup
awk '{print $1 "\t" $2}' "$POP_FILE_SOURCE" | sed 's/OUTGROUP/Outgroup/g' > SETS_pop.txt
echo "已生成 SETS_pop.txt (前5行):"
head -n 5 SETS_pop.txt

# 生成谱系水平配置 (SETS_lineage.txt): 第1列样本, 第3列谱系
awk '{print $1 "\t" $3}' "$POP_FILE_SOURCE" | sed 's/OUTGROUP/Outgroup/g' > SETS_lineage.txt
echo "已生成 SETS_lineage.txt (前5行):"
head -n 5 SETS_lineage.txt

# ==========================================
# 3. 创建进化树文件 (.nwk)
# ==========================================
echo "=== 3. 创建进化树文件 ==="

# 3.1 谱系树 (Lineage Tree)
# 结构: 外群 -> 西南云南 -> 西南贵州 -> (华中, 华东)
echo "(Outgroup,(Southwest-Yunnan,(Southwest-Guizhou,(Central,East))));" > tree_lineage.nwk
echo "已生成 tree_lineage.nwk"

# 3.2 精细群体树 (Population Tree)
# echo "已生成 tree_pop.nwk"

# ==========================================
# 4. 分析流程 A: 谱系水平 (Lineage Level)
# ==========================================
echo "=== 4. 开始运行谱系水平分析 (Lineage Level) ==="

# 4.1 运行 Dtrios
echo "运行 Dtrios..."
$DSUITE_BIN Dtrios \
    -t tree_lineage.nwk \
    -o result_lineage \
    "$INPUT_VCF" \
    SETS_lineage.txt

# 4.2 运行 Fbranch
echo "运行 Fbranch..."
$DSUITE_BIN Fbranch \
    tree_lineage.nwk \
    result_lineage_tree.txt \
    > result_lineage_fbranch.txt

# 4.3 绘图
echo "生成谱系热图..."
python "$DTOOLS_PY" \
    result_lineage_fbranch.txt \
    tree_lineage.nwk \
    --run-name result_lineage_plot

# ==========================================
# 5. 分析流程 B: 群体水平 (Population Level)
# ==========================================
echo "=== 5. 开始运行群体水平分析 (Population Level) ==="

# 5.1 运行 Dtrios
echo "运行 Dtrios (这可能需要几分钟)..."
$DSUITE_BIN Dtrios \
    -t tree_pop.nwk \
    -o result_pop \
    "$INPUT_VCF" \
    SETS_pop.txt

# 5.2 运行 Fbranch
echo "运行 Fbranch..."
$DSUITE_BIN Fbranch \
    tree_pop.nwk \
    result_pop_tree.txt \
    > result_pop_fbranch.txt

# 5.3 绘图
echo "生成群体热图..."
python "$DTOOLS_PY" \
    result_pop_fbranch.txt \
    tree_pop.nwk \
    --run-name result_pop_plot

# ==========================================
# 6. 完成
# ==========================================
echo "=========================================="
echo "所有分析已完成！"
echo "工作目录: $WORK_DIR"
echo "请下载查看以下图片文件:"
echo "1. 谱系分析图: result_lineage_plot.png (或 .svg)"
echo "2. 群体分析图: result_pop_plot.png (或 .svg)"
echo "=========================================="
