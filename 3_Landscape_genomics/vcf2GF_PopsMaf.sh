# ===========================
# Shell脚本部分：处理群体MAF
# ===========================

# --- 配置区域 ---
# 输入文件路径 (请仔细核对)
vcf_file="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.nomissing.recode.vcf.gz"
pop_file="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202samples.pop"

# 并发数控制 (防止电脑卡死，建议最多设置为 CPU 核心数的80%-90%)
MAX_JOBS=20

# 检查输入文件是否存在
if [ ! -f "$vcf_file" ]; then
    echo "Error: VCF file not found: $vcf_file"
    exit 1
fi
if [ ! -f "$pop_file" ]; then
    echo "Error: Population file not found: $pop_file"
    exit 1
fi

# 输出目录
output_dir="population_maf_results"
temp_dir="${output_dir}/temp_ids"
plink_binary_prefix="${output_dir}/master_binary"

# 确保输出和临时目录存在
mkdir -p "$output_dir"
mkdir -p "$temp_dir"

# ------------------------------------------------------------------
# Step 1: 预处理 VCF -> PLINK Binary
# ------------------------------------------------------------------
echo ">>> Step 1: 将 VCF 转换为 PLINK 二进制格式 (加速后续计算)..."
echo "    使用 --double-id 确保样本名同时作为 FID 和 IID，解决匹配问题。"

if [ ! -f "${plink_binary_prefix}.bed" ]; then
    plink --vcf "$vcf_file" \
          --make-bed \
          --allow-extra-chr \
          --keep-allele-order \
          --double-id \
          --out "$plink_binary_prefix" \
          --silent
    
    if [ $? -ne 0 ]; then
        echo "Error: VCF 转 PLINK 失败，请检查 VCF 格式。"
        exit 1
    fi
else
    echo "    检测到二进制文件已存在，跳过转换步骤。"
fi

# ------------------------------------------------------------------
# Step 2: 准备群体 ID 列表
# ------------------------------------------------------------------
echo ">>> Step 2: 提取并生成群体 ID 列表..."

# 提取唯一的群体名称
cut -f2 "$pop_file" | sort | uniq > "${temp_dir}/unique_groups.txt"

# 为每个群体生成 .ids 文件
# 格式要求：FamilyID IndividualID
# 由于 Step 1 用了 --double-id，所以这里两列都必须是 SampleID
while read -r group; do
    # 在 pop 文件中找该群体的样本，输出两列相同的 ID
    grep -w "$group" "$pop_file" | awk '{print $1, $1}' > "${temp_dir}/${group}.ids"
done < "${temp_dir}/unique_groups.txt"

# ------------------------------------------------------------------
# Step 3: 并行计算频率
# ------------------------------------------------------------------
echo ">>> Step 3: 开始并行计算各群体频率..."

# 定义处理函数
process_group() {
    local group=$1
    
    # 使用 --bfile 读取二进制文件，速度极快
    plink --bfile "$plink_binary_prefix" \
          --keep "${temp_dir}/${group}.ids" \
          --allow-extra-chr \
          --freq \
          --out "${output_dir}/${group}_population_maf" \
          --silent

    if [ $? -eq 0 ]; then
        echo "[完成] $group"
    else
        echo "[失败] $group (请检查日志 ${output_dir}/${group}_population_maf.log)"
    fi
}

# 循环提交任务
while read -r group; do
    process_group "$group" &
    
    # 并发控制：如果后台任务数 >= MAX_JOBS，则等待任意一个完成
    if [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; then
        wait -n
    fi
done < "${temp_dir}/unique_groups.txt"

# 等待所有剩余任务
wait
echo ">>> 所有群体频率计算完毕。"

# ====================================
# Python脚本：合并MAF结果 (高效版)
# ====================================
echo ">>> Step 4: 合并结果生成矩阵..."

python3 << 'EOF_PYTHON'
import pandas as pd
import os
import glob

def combine_maf_optimized(input_dir, output_file):
    print(f"正在读取 {input_dir} 下的频率文件...")
    
    # 获取所有 .frq 文件
    frq_files = glob.glob(os.path.join(input_dir, "*_population_maf.frq"))
    
    if not frq_files:
        print("Error: 未找到 .frq 文件，Step 3 可能全部失败。")
        return
    
    dfs = []
    
    for file_path in frq_files:
        # 从文件名提取群体名
        file_name = os.path.basename(file_path)
        pop_name = file_name.replace('_population_maf.frq', '')
        
        try:
            # 读取 PLINK 的 .frq 文件 (空格分隔)
            # 只需要 SNP 和 MAF 两列
            df = pd.read_csv(file_path, delim_whitespace=True, usecols=['SNP', 'MAF'])
            
            # 重命名 MAF 列为 群体名
            df.rename(columns={'MAF': pop_name}, inplace=True)
            
            # 设置 SNP 为索引，方便后续拼接
            df.set_index('SNP', inplace=True)
            
            dfs.append(df)
        except Exception as e:
            print(f"警告: 读取 {file_name} 失败 - {e}")

    if dfs:
        print("正在拼接矩阵 (这比循环 merge 快很多)...")
        # 横向拼接 (axis=1)
        combined_df = pd.concat(dfs, axis=1)
        
        # 填充 NaN (如果有群体缺失某些位点，填 0)
        combined_df.fillna(0, inplace=True)
        
        # 转置矩阵：变成 行=群体，列=SNP
        final_matrix = combined_df.T
        
        # 保存 CSV
        final_matrix.to_csv(output_file)
        print(f"成功！结果已保存到: {output_file}")
        print(f"矩阵维度: {final_matrix.shape}")
    else:
        print("Error: 没有有效数据被合并。")

if __name__ == "__main__":
    combine_maf_optimized("population_maf_results", "GF_PopsMaf.csv")
EOF_PYTHON

echo ">>> 脚本执行结束。"
