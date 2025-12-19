#!/bin/bash

# ===========================
# Shell脚本部分：处理群体MAF
# ===========================

# 输入文件和目录
vcf_file="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.nomissing.recode.vcf.gz"
pop_file="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202samples.pop"

# 输出目录
output_dir="population_maf_results"
temp_dir="${output_dir}/temp_ids"
plink_binary_prefix="${output_dir}/master_binary" # 中间二进制文件前缀

# 并发控制：设置最大同时运行的任务数（根据你的CPU核心数调整，例如 10）
MAX_JOBS=20

# 确保输出和临时目录存在
mkdir -p "$output_dir"
mkdir -p "$temp_dir"

echo "Step 1: 预处理 - 将 VCF 转换为 PLINK 二进制格式以加速读取..."
# 这一步只做一次，后续读取速度提升 10 倍以上
if [ ! -f "${plink_binary_prefix}.bed" ]; then
    plink --vcf "$vcf_file" \
          --make-bed \
          --allow-extra-chr \
          --keep-allele-order \
          --const-fid \
          --out "$plink_binary_prefix" \
          --silent
fi

echo "Step 2: 提取群体信息..."
cut -f2 "$pop_file" | sort | uniq > "${temp_dir}/unique_groups.txt"

# 为每个独特的群体生成包含个体ID的文件
while read -r group; do
    # 假设 pop 文件格式为: SampleID PopID
    # Plink keep 文件需要两列: FamilyID IndividualID
    # 如果你的 VCF 没有 FamilyID，Plink 通常将 FID 视为 IID 或 0。
    # 这里我们生成两列相同的 ID 以防万一 (FID IID)
    awk -v g="$group" '$2 == g {print $1, $1}' "$pop_file" > "${temp_dir}/${group}.ids"
done < "${temp_dir}/unique_groups.txt"

echo "Step 3: 开始并行计算频率..."

# 定义计算函数
process_group() {
    local group=$1
    # 直接基于转换好的二进制文件计算频率，无需生成中间bed文件
    plink --bfile "$plink_binary_prefix" \
          --keep "${temp_dir}/${group}.ids" \
          --allow-extra-chr \
          --freq \
          --out "${output_dir}/${group}_population_maf" \
          --silent

    if [ $? -eq 0 ]; then
        echo "[Success] $group finished."
    else
        echo "[Error] $group failed."
    fi
}

# 循环并控制并发
while read -r group; do
    process_group "$group" &
    
    # 任务计数控制
    if [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; then
        wait -n
    fi
done < "${temp_dir}/unique_groups.txt"

# 等待剩余任务完成
wait

# 清理中间二进制大文件（可选，如果不需要保留）
# rm "${plink_binary_prefix}.*"

echo "所有群体频率计算完毕。"

# ====================================
# 嵌入第二个Python脚本：合并MAF结果 (优化版)
# ====================================
# 注：省略了第一个简单的分析脚本，直接进行核心的合并步骤

python3 << 'EOF2'
import pandas as pd
import os
import glob

def combine_maf_optimized(input_dir, output_file):
    print("===== 开始合并MAF结果 (优化版) =====")
    
    # 获取所有 .frq 文件
    frq_files = glob.glob(os.path.join(input_dir, "*_population_maf.frq"))
    
    if not frq_files:
        print("未找到任何 .frq 文件用于合并。")
        return
    
    print(f"找到 {len(frq_files)} 个群体文件，开始读取...")
    
    dfs = []
    
    for file_path in frq_files:
        # 提取群体名称
        file_name = os.path.basename(file_path)
        population_name = file_name.replace('_population_maf.frq', '')
        
        try:
            # Plink .frq 文件通常是空格分隔
            # 列名: CHR SNP A1 A2 MAF NCHROBS
            # 我们只需要 SNP 和 MAF
            df = pd.read_csv(file_path, sep=r'\s+', usecols=['SNP', 'MAF'], engine='c')
            
            # 将 SNP 设为索引，列名改为群体名
            df.set_index('SNP', inplace=True)
            df.rename(columns={'MAF': population_name}, inplace=True)
            
            dfs.append(df)
            
        except Exception as e:
            print(f"读取文件 {file_name} 失败: {e}")

    if dfs:
        print("正在拼接矩阵 (这可能需要一点时间)...")
        # 使用 concat 一次性拼接，比循环 merge 快得多
        # axis=1 表示横向拼接 (按列)
        combined_maf = pd.concat(dfs, axis=1)
        
        # 如果需要转置 (行=群体, 列=SNP)
        print("正在转置并保存...")
        combined_maf_T = combined_maf.T
        
        # 保存
        combined_maf_T.to_csv(output_file)
        print(f"===== 成功！结果已保存到 {output_file} =====")
        print(f"矩阵维度: {combined_maf_T.shape}")
    else:
        print("没有有效的数据被合并。")

def main():
    input_directory = "population_maf_results"
    output_csv = "GF_PopsMaf.csv"
    combine_maf_optimized(input_directory, output_csv)

if __name__ == "__main__":
    main()
EOF2

echo "脚本执行结束。"
