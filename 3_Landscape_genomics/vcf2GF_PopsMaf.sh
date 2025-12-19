# ===========================
# Shell脚本部分：处理群体MAF
# ===========================

# 输入文件路径
vcf_file="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.nomissing.recode.vcf.gz"
pop_file="/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202samples.pop"

# 检查输入文件
if [ ! -f "$vcf_file" ] || [ ! -f "$pop_file" ]; then
    echo "【错误】找不到输入文件，请检查路径！"
    exit 1
fi

# 输出目录
output_dir="population_maf_results"
temp_dir="${output_dir}/temp_ids"
plink_binary_prefix="${output_dir}/master_binary" 
MAX_JOBS=20

mkdir -p "$output_dir"
mkdir -p "$temp_dir"

echo "Step 1: 预处理 - 将 VCF 转为 PLINK 格式 (使用 double-id 修复ID匹配)..."

# 【关键修改】使用 --double-id，确保 FID 和 IID 都是样本名
plink --vcf "$vcf_file" \
      --make-bed \
      --allow-extra-chr \
      --keep-allele-order \
      --double-id \
      --out "$plink_binary_prefix" \
      --silent

# 检查是否成功
if [ ! -f "${plink_binary_prefix}.fam" ]; then
    echo "【错误】PLINK 转换失败，未生成 .fam 文件。"
    exit 1
fi

echo "Step 2: 提取群体信息..."
cut -f2 "$pop_file" | sort | uniq > "${temp_dir}/unique_groups.txt"

# 生成 ID 列表
while read -r group; do
    # 生成 FID IID 两列，均使用样本名
    awk -v g="$group" '$2 == g {print $1, $1}' "$pop_file" > "${temp_dir}/${group}.ids"
done < "${temp_dir}/unique_groups.txt"

# ================= Debug 环节 =================
# 在大规模运行前，先进行一次 ID 匹配测试，防止刷屏报错
echo "Running ID match check..."
first_group=$(head -n 1 "${temp_dir}/unique_groups.txt")
head_fam=$(head -n 3 "${plink_binary_prefix}.fam" | awk '{print $1, $2}')
head_ids=$(head -n 3 "${temp_dir}/${first_group}.ids")

echo "--- Debug Info ---"
echo "PLINK 数据中的 ID (前3行 FID IID):"
echo "$head_fam"
echo "我们生成的筛选 ID (前3行 FID IID):"
echo "$head_ids"
echo "------------------"
# ==============================================

echo "Step 3: 开始并行计算频率..."

process_group() {
    local group=$1
    plink --bfile "$plink_binary_prefix" \
          --keep "${temp_dir}/${group}.ids" \
          --allow-extra-chr \
          --freq \
          --out "${output_dir}/${group}_population_maf" \
          --silent

    if [ $? -eq 0 ]; then
        echo "[Success] $group finished."
    else
        echo "[Error] $group failed (Check logs)."
    fi
}

while read -r group; do
    process_group "$group" &
    if [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; then
        wait -n
    fi
done < "${temp_dir}/unique_groups.txt"

wait
echo "所有群体频率计算完毕。"

# ====================================
# Python脚本：合并MAF结果
# ====================================

python3 << 'EOF2'
import pandas as pd
import os
import glob

def combine_maf_optimized(input_dir, output_file):
    print("===== 开始合并MAF结果 =====")
    frq_files = glob.glob(os.path.join(input_dir, "*_population_maf.frq"))
    
    if not frq_files:
        print("错误：未找到任何 .frq 文件。")
        return
    
    print(f"找到 {len(frq_files)} 个文件，开始合并...")
    dfs = []
    
    for file_path in frq_files:
        file_name = os.path.basename(file_path)
        # 移除后缀得到群体名
        population_name = file_name.replace('_population_maf.frq', '')
        try:
            df = pd.read_csv(file_path, sep=r'\s+', usecols=['SNP', 'MAF'], engine='c')
            df.set_index('SNP', inplace=True)
            df.rename(columns={'MAF': population_name}, inplace=True)
            dfs.append(df)
        except Exception as e:
            print(f"读取失败: {file_name}")

    if dfs:
        combined_maf = pd.concat(dfs, axis=1)
        combined_maf.T.to_csv(output_file)
        print(f"成功！结果已保存到 {output_file}")
    else:
        print("无数据合并。")

if __name__ == "__main__":
    combine_maf_optimized("population_maf_results", "GF_PopsMaf.csv")
EOF2

echo "脚本结束。"
