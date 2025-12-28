#!/bin/bash

# ================= 1. 全局配置区域 =================
JAR_PATH="/home/vensin/software/FitCoal1.3/FitCoal.jar"
MU_PER_KB="0.00000419"
GENOME_LEN_KB="693703"
GEN_TIME="6"

# 并行控制
MAX_JOBS=24       # 同时运行的最大任务数
BOOT_REPEATS=50   # 生成多少个 Bootstrap 样本 ( 50-100 次用于画置信区间)
FITCOAL_REPEATS=20 # FitCoal 计算时的重复次数 (Bootstrap 仅需看分布，20次足够，追求速度)

# 目录设置
INPUT_DIR="bootstrap/bootstrap_inputs"   # 存放生成的 SFS 文件
OUTPUT_DIR="bootstrap/bootstrap_results" # 存放计算结果
TABLES_DIR="tables"            # 原始 tables 目录位置

# ================= 2. 环境准备 =================
echo ">>> [Step 1] 环境初始化..."

# 创建必要的目录
if [ ! -d "$INPUT_DIR" ]; then mkdir -p "$INPUT_DIR"; fi
if [ ! -d "$OUTPUT_DIR" ]; then mkdir -p "$OUTPUT_DIR"; fi

# 链接 Tables 目录 (避免重复造表)
if [ ! -d "$OUTPUT_DIR/tables" ]; then
    if [ -d "$TABLES_DIR" ]; then
        echo "    -> 建立 tables 软链接..."
        ln -s "$(realpath $TABLES_DIR)" "$OUTPUT_DIR/tables"
    else
        echo "    -> 警告: 未找到原始 tables 目录，程序将自动重新生成。"
        mkdir -p "$OUTPUT_DIR/tables"
    fi
fi

# ================= 3. 嵌入 Python 脚本生成数据 =================
echo ">>> [Step 2] 生成 Bootstrap SFS 数据..."

# 使用 Here Document 生成临时的 python 脚本
cat << 'EOF' > generate_boot_temp.py
import numpy as np
import glob
import os
import sys

# 接收外部参数
input_suffix = ".fitcoal.txt"
n_bootstraps = int(sys.argv[1])
output_dir = sys.argv[2]

sfs_files = glob.glob(f"*{input_suffix}")
print(f"    -> 发现原始文件: {len(sfs_files)} 个")

for sfs_file in sfs_files:
    pop_name = sfs_file.replace(input_suffix, "")
    
    # 读取原始 SFS
    with open(sfs_file, 'r') as f:
        content = f.read().strip()
        numbers = [int(x) for x in content.split() if x.isdigit()]
        original_sfs = np.array(numbers)
    
    total_snps = np.sum(original_sfs)
    if total_snps == 0: continue
    probs = original_sfs / total_snps
    
    # 生成 Bootstrap
    for i in range(1, n_bootstraps + 1):
        boot_sfs = np.random.multinomial(total_snps, probs)
        out_filename = os.path.join(output_dir, f"{pop_name}_boot_{i}{input_suffix}")
        with open(out_filename, 'w') as out_f:
            out_f.write("\t".join(map(str, boot_sfs)) + "\n")
EOF

# 运行 Python 生成数据
python3 generate_boot_temp.py "$BOOT_REPEATS" "$INPUT_DIR"
# 删除临时 Python 脚本
rm generate_boot_temp.py
echo "    -> Bootstrap 数据生成完毕！存放于 $INPUT_DIR"

# ================= 4. 智能排序与并行计算 =================
echo ">>> [Step 3] 开始并行计算 (优先处理大样本群体)..."

# --- 核心逻辑：按样本量排序 ---
# 1. awk 计算每个文件的 SFS 总和 (样本量)
# 2. sort -rn 按数字倒序排列 (大的在前)
# 3. cut 提取文件名
# 这样 样本量大的任务会排在最前面，样本量小的任务排在后面
echo "    -> 正在计算任务优先级..."
SORTED_FILES=$(awk '{sum=0; for(i=1; i<=NF; i++) sum+=$i; print sum, FILENAME}' $INPUT_DIR/*.fitcoal.txt | sort -rn | cut -d' ' -f2)

count=0
total_files=$(ls $INPUT_DIR/*.fitcoal.txt | wc -l)

for input_path in $SORTED_FILES; do
    ((count++))
    filename=$(basename "$input_path")
    task_name="${filename%.fitcoal.txt}"
    
    # 检查是否已运行 (断点续传)
    if [ -f "$OUTPUT_DIR/${task_name}.fitcoal.output" ]; then
        echo "    [$count/$total_files] 跳过已完成: $task_name"
        continue
    fi
    
    # 任务控制：后台任务数 >= MAX_JOBS 时等待
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 1
    done
    
    echo "    [$count/$total_files] 启动: $task_name (后台)"
    
    # --- 运行 FitCoal ---
    # output 指定到子目录
    # log 也重定向到子目录
    nohup java -Xmx2g -cp "$JAR_PATH" FitCoal.calculate.SinglePopDecoder \
        -table "$OUTPUT_DIR/tables/" \
        -input "$input_path" \
        -output "$OUTPUT_DIR/${task_name}.fitcoal.output" \
        -mutationRate "$MU_PER_KB" \
        -generationTime "$GEN_TIME" \
        -genomeLength "$GENOME_LEN_KB" \
        -repeats "$FITCOAL_REPEATS" > "$OUTPUT_DIR/${task_name}.log" 2>&1 &
        
done

echo ""
echo ">>> 所有任务已提交！正在等待后台计算结束..."
wait
echo ">>> ========================================="
echo ">>> 全部完成！"
echo ">>> 结果保存在: $OUTPUT_DIR"
echo ">>> ========================================="
