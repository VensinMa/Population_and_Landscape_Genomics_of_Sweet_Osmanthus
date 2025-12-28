cat > run_fitcoal_parallel.sh << 'EOF'
#!/bin/bash

# ================= 配置区域 =================
JAR_PATH="/home/vensin/software/FitCoal1.3/FitCoal.jar"

# 参数设置
MU_PER_KB="0.00000419"
GENOME_LEN_KB="693703"
GEN_TIME="6"

# ================= 开始计算 =================
echo ">>> 开始启动并行计算..."
echo ">>> 注意：所有群体/谱系将同时计算，日志输出到同名 .log 文件中。"
echo ">>> 脚本将保持运行状态，直到所有任务结束 (Ctrl+C 同时终止所有任务)。"
echo ""

# 确保 tables 存在
if [ ! -d "tables" ]; then
    mkdir tables
fi

# 循环启动任务
for input_file in *.fitcoal.txt; do
    pop_name="${input_file%.fitcoal.txt}"
    
    echo "[启动任务] 群体: $pop_name"
    
    # 核心命令解释：
    java -Xmx10g -cp "$JAR_PATH" FitCoal.calculate.SinglePopDecoder \
        -table tables/ \
        -input "$input_file" \
        -output "${pop_name}.fitcoal.output" \
        -mutationRate "$MU_PER_KB" \
        -generationTime "$GEN_TIME" \
        -genomeLength "$GENOME_LEN_KB" \
        -repeats 100 > "${pop_name}.log" 2>&1 &
        
    # 记录进程ID (可选，方便看)
    pid=$!
    echo "   -> PID: $pid (日志: ${pop_name}.log)"
done

echo ""
echo ">>> 所有任务已提交，正在并行计算中... (请勿关闭终端)"
echo ">>> 你可以另开一个终端用 'tail -f *.log' 查看进度"

# wait 会卡住脚本，直到上面所有用 & 启动的子进程都结束
wait

echo ""
echo "========================================"
echo "所有群体计算全部完成！"
echo "========================================"
EOF
