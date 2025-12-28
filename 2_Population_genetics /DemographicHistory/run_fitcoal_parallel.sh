cat > run_fitcoal_parallel.sh << 'EOF'
#!/bin/bash

# ================= 配置区域 =================
# FitCoal.jar 的路径 
JAR_PATH="/home/vensin/software/FitCoal1.3/FitCoal.jar"

# 参数设置
MU_PER_KB="0.00000419"
GENOME_LEN_KB="693703"
GEN_TIME="6"

# ================= 开始计算 =================
echo "开始并行计算..."

# 循环处理目录下的所有 .fitcoal.txt 文件
for input_file in *.fitcoal.txt; do
    # 获取群体名
    pop_name="${input_file%.fitcoal.txt}"
    
    echo "----------------------------------------"
    echo "正在后台启动群体: $pop_name"
    echo "----------------------------------------"
    
    # 运行 FitCoal
    # 确保 tables 文件夹存在
    mkdir tables
    nohup java -Xmx10g -cp "$JAR_PATH" FitCoal.calculate.SinglePopDecoder \
        -table tables/ \
        -input "$input_file" \
        -output "${pop_name}.fitcoal.output" \
        -mutationRate "$MU_PER_KB" \
        -generationTime "$GEN_TIME" \
        -genomeLength "$GENOME_LEN_KB" \
        -repeats 100 > "${pop_name}.log" 2>&1 &
        
done

# 等待所有后台任务完成
echo "所有任务已提交，正在计算中..."
wait
echo "========================================"
echo "所有群体计算全部完成！"
echo "========================================"
EOF
