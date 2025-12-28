cat > run_fitcoal.sh << 'EOF'
#!/bin/bash

# ================= 配置区域 =================
# FitCoal.jar 的路径 (根据你之前的 ll 信息)
JAR_PATH="/home/vensin/software/FitCoal1.3/FitCoal.jar"

# 参数设置 (FitCoal输入数据需要换算)
# 4.19e-9 * 1000
MU_PER_KB="0.00000419"
# 693703244 / 1000
GENOME_LEN_KB="693703"
GEN_TIME="6"

# ================= 开始计算 =================
# 循环处理目录下的所有 .sfs.txt 文件
for input_file in *.sfs.txt; do
    # 获取群体名 (如 Central)
    pop_name="${input_file%.sfs.txt}"
    
    echo "----------------------------------------"
    echo "正在运行群体: $pop_name"
    echo "突变率: $MU_PER_KB (per KB)"
    echo "----------------------------------------"
    
    # 运行 FitCoal
    # -table tables/ : 必须指定，程序会自动创建
    # -repeats 100    : 重复计算 100 次取最优
    java -Xmx40g -cp "$JAR_PATH" FitCoal.calculate.SinglePopDecoder \
        -table tables/ \
        -input "$input_file" \
        -output "${pop_name}.fitcoal.output" \
        -mutationRate "$MU_PER_KB" \
        -generationTime "$GEN_TIME" \
        -genomeLength "$GENOME_LEN_KB" \
        -repeats 100
        
    echo "群体 $pop_name 计算完成！"
done
EOF
