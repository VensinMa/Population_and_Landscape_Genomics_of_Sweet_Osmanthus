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

if [ ! -f "$POP_FILE_SOURCE" ]; then
    echo "错误：找不到群体映射文件 $POP_FILE_SOURCE"
    exit 1
fi

# ==========================================
# 2. 准备所有配置文件 (串行执行，确保文件就绪)
# ==========================================
echo "=== 2. 生成所有配置文件 ==="

# 生成 SETS 文件
awk '{print $1 "\t" $2}' "$POP_FILE_SOURCE" | sed 's/OUTGROUP/Outgroup/g' > SETS_pop.txt
awk '{print $1 "\t" $3}' "$POP_FILE_SOURCE" | sed 's/OUTGROUP/Outgroup/g' > SETS_lineage.txt

# 生成树文件
echo "(Outgroup,(Southwest-Yunnan,(Southwest-Guizhou,(Central,East))));" > tree_lineage.nwk
# 注意：这里请填入您之前确认的精细群体树结构
echo "(Outgroup:0.0753917,((LCJ:0.08767953,DRS:0.08844153):0.0311433,((JMX:0.141584,(ZJS:0.14525122,RX:0.14878245):0.00067956):0.0109852,(HYX:0.169961,(HJLL:0.174504357,((ZJP:0.18144733,(EJ:0.180232354,(DA:0.18386679,DST:0.177460344):0.00867956):0.00221941):0.000558854,((XC:0.176513,YX:0.1658082):0.0089426,(SLZ:0.1787979,((SXK:0.1661001,SFZ:0.16691828):0.01397163,(SK:0.178496493,((XNF:0.17725563,(LQ:0.17747174,CP:0.178287382):0.000633262):6.91048e-05,((JD:0.17721253,QDH:0.17468136):0.000682841,((CX:0.17569913,(LX:0.1746603,(SL:0.166422,WYL:0.16978951):0.0078507):0.00152176):0.000461085,(LS:0.17900896,(YK:0.1763265,(YZY:0.17642056,(GHX:0.17999454,(FZA:0.17531431,ST:0.2278669):0.0051006):0.00136227):0.00117963):0.000680322):5.54522e-05):5.7507e-05):0.000662918):0.000440513):0.00112392):0.000339692):0.00214137):0.00408189):0.00752432):0.0105002):0.03444564):0.0405461):0.0471352);" > tree_pop.nwk

echo "配置文件准备完毕，开始 Dsuite 分析..."
echo "----------------------------------------"

# ==========================================
# 3. 定义函数：谱系水平分析
# ==========================================
run_lineage_analysis() {
    echo "[谱系分析] 开始运行 Dtrios..."
    $DSUITE_BIN Dtrios -t tree_lineage.nwk -o result_lineage "$INPUT_VCF" SETS_lineage.txt > lineage_dtrios.log 2>&1
    
    echo "[谱系分析] Dtrios 完成，开始 Fbranch..."
    $DSUITE_BIN Fbranch tree_lineage.nwk result_lineage_tree.txt > result_lineage_fbranch.txt
    
    echo "[谱系分析] Fbranch 完成，开始绘图..."
    python "$DTOOLS_PY" result_lineage_fbranch.txt tree_lineage.nwk --run-name result_lineage_plot > lineage_plot.log 2>&1
    
    echo "[谱系分析] 全部完成！"
}

# ==========================================
# 4. 定义函数：群体水平分析
# ==========================================
run_population_analysis() {
    echo "[群体分析] 开始运行 Dtrios (耗时较长)..."
    $DSUITE_BIN Dtrios -t tree_pop.nwk -o result_pop "$INPUT_VCF" SETS_pop.txt > pop_dtrios.log 2>&1
    
    echo "[群体分析] Dtrios 完成，开始 Fbranch..."
    $DSUITE_BIN Fbranch tree_pop.nwk result_pop_tree.txt > result_pop_fbranch.txt
    
    echo "[群体分析] Fbranch 完成，开始绘图..."
    python "$DTOOLS_PY" result_pop_fbranch.txt tree_pop.nwk --run-name result_pop_plot > pop_plot.log 2>&1
    
    echo "[群体分析] 全部完成！"
}

# ==========================================
# 5. 并行执行任务
# ==========================================

# 在后台启动任务 A
run_lineage_analysis &
PID_LINEAGE=$!
echo "谱系任务已在后台启动 (PID: $PID_LINEAGE)"

# 在后台启动任务 B
run_population_analysis &
PID_POP=$!
echo "群体任务已在后台启动 (PID: $PID_POP)"

# ==========================================
# 6. 等待所有任务结束
# ==========================================
echo "正在等待两个任务完成..."
wait $PID_LINEAGE
EXIT_CODE_LINEAGE=$?

wait $PID_POP
EXIT_CODE_POP=$?

# ==========================================
# 7. 检查结果
# ==========================================
echo "----------------------------------------"
if [ $EXIT_CODE_LINEAGE -eq 0 ]; then
    echo "✅ 谱系分析成功完成。"
else
    echo "❌ 谱系分析失败，请查看 lineage_dtrios.log 或 lineage_plot.log"
fi

if [ $EXIT_CODE_POP -eq 0 ]; then
    echo "✅ 群体分析成功完成。"
else
    echo "❌ 群体分析失败，请查看 pop_dtrios.log 或 pop_plot.log"
fi

echo "=========================================="
echo "工作目录: $WORK_DIR"
echo "请查看生成的 result_lineage_plot*.png 和 result_pop_plot*.png"
