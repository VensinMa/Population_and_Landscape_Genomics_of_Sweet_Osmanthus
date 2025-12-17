#  bed格式
## plink --vcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.recode.vcf.gz \
##     --make-bed   --out /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned  --keep-allele-order  --allow-extra-chr

# 1. 创建输出目录
mkdir -p /home/vensin/workspace/snpcalling_wild/12.population_genetics/pca
cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/pca

# 2. 运行 GCTA 计算 GRM
# --bfile: 指定输入文件前缀 (必须是 PLINK 的 .bed/.bim/.fam 格式)
# --make-grm: 构建矩阵
# --autosome: 仅使用常染色体 (自动排除性染色体)
# --thread-num 10: 建议加上多线程，加快计算速度
# --out: 输出前缀 (这里命名为 202_samples_pca)

gcta64 --bfile /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned \
       --make-grm \
       --thread-num 8 \
       --out 202_samples_gcta_pca

# 3. 运行 PCA (提取前 3 个主成分)
# --grm: 输入上一步生成的 GRM 前缀
# --pca 3: 提取前 3 个主成分
# --out: 输出文件名

gcta64 --grm 202_samples_gcta_pca \
       --pca 3 \
       --out 202_samples_gcta_pca_result
