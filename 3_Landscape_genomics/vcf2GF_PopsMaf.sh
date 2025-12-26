#!/bin/bash

# ===========================
# 优化版 Shell脚本 v2：修复 ID 匹配问题
# ===========================

# 输入文件和目录
vcf_file="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.nomissing.recode.vcf.gz"
pop_file="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202samples.pop"

# 输出目录
output_dir="population_maf_results"
temp_dir="${output_dir}/temp_ids"
global_bed="${output_dir}/global_dataset" # 全局二进制文件前缀

# 并发控制：同时运行的最大任务数
MAX_JOBS=10

mkdir -p "$output_dir"
mkdir -p "$temp_dir"

echo "=== 第一步：预处理 ==="
# 1. 提取群体信息
cut -f2 "$pop_file" | sort | uniq > "${temp_dir}/unique_groups.txt"

# 2. 生成ID列表
# 格式：SampleID SampleID (为了匹配 --double-id)
while read -r group; do
    awk -v g="$group" '$2 == g {print $1, $1}' "$pop_file" > "${temp_dir}/${group}.ids"
done < "${temp_dir}/unique_groups.txt"

# 3. [关键修正] 将 VCF 转为全局 PLINK 二进制文件
# 增加了 --double-id 参数，确保 .fam 文件中的 FID 和 IID 都是样本名
if [ ! -f "${global_bed}.bed" ]; then
    echo "正在将 VCF 转换为全局 PLINK 二进制文件 (包含 --double-id)..."
    plink --vcf "$vcf_file" \
          --make-bed \
          --double-id \
          --allow-extra-chr \
          --keep-allele-order \
          --out "$global_bed" \
          --set-missing-var-ids @:# \
          --silent
    
    # [Debug] 检查一下生成的 ID 格式，输出前5行给用户看
    echo "--- 检查生成的 .fam 文件前5行 ---"
    head -n 5 "${global_bed}.fam"
    echo "--------------------------------"
else
    echo "全局 PLINK 文件已存在，跳过转换。"
fi

echo "=== 第二步：并行计算各群体 MAF ==="

# 任务计数器
count=0

while read -r group; do
    (
        # 直接读取全局 bed 文件，只计算频率
        plink --bfile "$global_bed" \
              --keep "${temp_dir}/${group}.ids" \
              --freq \
              --allow-extra-chr \
              --out "${output_dir}/${group}_population_maf" \
              --silent

        if [ -f "${output_dir}/${group}_population_maf.frq" ]; then
            echo "[成功] 群体 $group 处理完毕"
        else
            # 如果失败，打印一下对应的 ID 文件内容，方便调试
            echo "[失败] 群体 $group 计算出错。请检查 ${temp_dir}/${group}.ids 是否与 .fam 文件匹配"
        fi
    ) &

    # 并发控制逻辑
    ((count++))
    if [ $((count % MAX_JOBS)) -eq 0 ]; then
        wait # 每启动 MAX_JOBS 个任务后等待它们完成
    fi

done < "${temp_dir}/unique_groups.txt"

wait
echo "所有群体 MAF 计算完毕。"

# ====================================
# 嵌入优化后的 Python脚本：合并结果
# ====================================

python3 << 'EOF'
import pandas as pd
import os
import glob

def efficient_combine_maf(input_dir, output_file):
    print("\n===== 开始合并 MAF 结果 (Pandas 优化版) =====")
    
    frq_files = glob.glob(os.path.join(input_dir, "*_population_maf.frq"))
    
    if not frq_files:
        print("错误：未找到 .frq 文件。这说明上一步 PLINK 计算全部失败。")
        return

    data_frames = []
    
    for file_path in frq_files:
        # 获取群体名
        file_name = os.path.basename(file_path)
        population_name = file_name.replace('_population_maf.frq', '')
        
        try:
            # 读取 SNP 和 MAF 列
            df = pd.read_csv(file_path, sep=r'\s+', usecols=['SNP', 'MAF'], engine='c')
            
            # 设置 SNP 为索引
            df.set_index('SNP', inplace=True)
            df.rename(columns={'MAF': population_name}, inplace=True)
            
            data_frames.append(df)
            
        except Exception as e:
            print(f"读取文件 {file_name} 失败: {e}")

    if not data_frames:
        print("没有有效数据被读取。")
        return

    print(f"正在合并 {len(data_frames)} 个群体的数据...")
    
    # 使用 concat 横向拼接
    combined_maf = pd.concat(data_frames, axis=1, join='outer')
    
    # 填充缺失值 (可选，这里填0表示该群体该位点无多态或缺失)
    # combined_maf.fillna(0, inplace=True) 

    # 转置
    final_df = combined_maf.T
    
    # 保存
    final_df.to_csv(output_file)
    print(f"===== 合并完成 =====")
    print(f"输出文件: {output_file}")
    print(f"矩阵维度: {final_df.shape} (行:群体, 列:SNPs)")

if __name__ == "__main__":
    efficient_combine_maf("population_maf_results", "GF_PopsMaf.csv")
EOF

echo "脚本执行结束。"
