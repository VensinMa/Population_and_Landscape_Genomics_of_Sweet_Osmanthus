#!/bin/bash
# vcf2dis_bootstrap_nj_tree_Parallel_v2.sh - 修复函数传递问题的版本
# 用法: ./vcf2dis_bootstrap_nj_tree_Parallel_v2.sh <input.vcf> <bootstrap_replicates> <threads>

# ==============================================================================
# 软件路径配置
# ==============================================================================
VCF2DIS_BIN="/home/vensin/software/VCF2Dis-1.55/bin/VCF2Dis"
FNEIGHBOR_BIN="/home/vensin/anaconda3/envs/embassy-phylip/bin/fneighbor"
FCONSENSE_BIN="/home/vensin/anaconda3/envs/embassy-phylip/bin/fconsense"
ADD_BOOTSTRAP_PL="/home/vensin/software/VCF2Dis-1.55/bin/percentageboostrapTree.pl"
# ==============================================================================

# 检查参数
if [ "$#" -ne 3 ]; then
    echo "错误: 需要三个参数"
    echo "用法: $0 <input.vcf.gz> <bootstrap_times> <threads>"
    exit 1
fi

INPUT_VCF=$1
NN=$2
NUM_CORES=$3

# 检查输入文件
if [ ! -f "$INPUT_VCF" ]; then
    echo "错误: 输入文件不存在: $INPUT_VCF"
    exit 1
fi

# 检查软件
for tool in "$VCF2DIS_BIN" "$FNEIGHBOR_BIN" "$FCONSENSE_BIN" "$ADD_BOOTSTRAP_PL"; do
    if [ ! -f "$tool" ]; then
        echo "错误: 找不到软件: $tool"
        exit 1
    fi
done

echo "=== 开始 Bootstrap NJ 树分析 ==="
echo "输入: $INPUT_VCF"
echo "次数: $NN"
echo "线程: $NUM_CORES"

# 创建临时目录
TMP_DIR="bootstrap_tmp_$(date +%s)"
mkdir -p "$TMP_DIR"

# 定义核心函数
run_bootstrap() {
    # 在函数内部接收参数
    local id=$1
    local tmp_dir=$2
    local vcf=$3
    # 在函数内部接收软件路径 (通过环境变量传入)
    local vcf2dis=$4
    local fneighbor=$5

    # 1. 生成距离矩阵
    "$vcf2dis" -InPut "$vcf" -OutPut "$tmp_dir/p_dis_${id}.mat" -Rand 0.25 > /dev/null 2>&1
    
    if [ ! -f "$tmp_dir/p_dis_${id}.mat" ]; then
        echo "Error: Matrix ${id} failed"
        return 1
    fi
    
    # 2. 构建 NJ 树
    "$fneighbor" -datafile "$tmp_dir/p_dis_${id}.mat" \
                 -outfile "$tmp_dir/tree.info_${id}.txt" \
                 -matrixtype s \
                 -treetype n \
                 -outtreefile "$tmp_dir/tree_${id}.tre" > /dev/null 2>&1
}

# === 关键修改 ===
# 1. 获取函数定义的字符串
FUNC_BODY=$(declare -f run_bootstrap)

# 2. 导出所有必要的软件路径变量，以便 parallel 子shell 可以看到
export VCF2DIS_BIN FNEIGHBOR_BIN

echo "步骤1/3: 并行生成 Bootstrap 树..."

# 3. 执行 parallel
# 技巧: 我们把函数定义 ($FUNC_BODY) 直接放在命令序列的最前面
# 这样子 shell 一启动就会先"学会"这个函数，然后再执行它
seq 1 "$NN" | parallel -j "$NUM_CORES" --progress --bar \
    "$FUNC_BODY; run_bootstrap {} $TMP_DIR $INPUT_VCF $VCF2DIS_BIN $FNEIGHBOR_BIN"

# 检查结果
SUCCESS_COUNT=$(ls "$TMP_DIR"/tree_*.tre 2>/dev/null | wc -l)
echo ""
echo "成功生成: $SUCCESS_COUNT / $NN"

if [ "$SUCCESS_COUNT" -eq 0 ]; then
    echo "严重错误: 没有任何树生成，请检查上面的报错信息。"
    rm -rf "$TMP_DIR"
    exit 1
fi

# 合并树
echo "步骤2/3: 合并树文件..."
cat "$TMP_DIR"/tree_*.tre > "$TMP_DIR/alltree_merge.tre"

# 生成一致树
echo "步骤3/3: 生成一致树 (Consensus Tree)..."
"$FCONSENSE_BIN" -intreefile "$TMP_DIR/alltree_merge.tre" \
          -outfile "$TMP_DIR/consensus.info" \
          -outtreefile "$TMP_DIR/consensus.tre" \
          -treeprint Y > /dev/null 2>&1

# 添加支持率
echo "添加 Bootstrap 支持率到节点..."
perl "$ADD_BOOTSTRAP_PL" "$TMP_DIR/alltree_merge.tre" "$SUCCESS_COUNT" Final_boostrap.tre

# 复制结果
cp "$TMP_DIR/consensus.tre" .
cp "$TMP_DIR/consensus.info" .

echo "=== 分析完成 ==="
echo "结果文件: consensus.tre"
