#####################################  SNP ##############################################################
# 1. 创建并进入新目录
mkdir -p mkdir -p /home/vensin/workspace/snpcalling_wild/12.Vcftools_filter_no_outgroup/snp
cd /home/vensin/workspace/snpcalling_wild/12.Vcftools_filter_no_outgroup/snp

# 2. 运行 VCFtools (SNP)
# 新增: --remove-indv 参数重复 6 次以排除特定个体
# 注意: 输出文件名改为 202_samples...
vcftools --gzvcf /home/data/10.gatk_variantfiltration/SNP/snp_filtered.vcf.gz \
  --remove-indv O_DSMX \
  --remove-indv O_MXL \
  --remove-indv O_MZGH \
  --remove-indv O_WMMX \
  --remove-indv O_XYYG \
  --remove-indv O_ZNMX \
  --min-alleles 2 --max-alleles 2 \
  --minGQ 10 --minQ 30 --min-meanDP 6 \
  --max-missing 0.9 --maf 0.05 \
  --recode --recode-INFO-all \
  --stdout \
  2> 202_samples_snp.log \
  | bgzip -@ 8 > 202_samples_snp_filtered.recode.vcf.gz

# 3. 建立索引
tabix -p vcf 202_samples_snp_filtered.recode.vcf.gz

#####################################  INDEL ##############################################################
# 1. 创建并进入新目录
mkdir -p /home/vensin/workspace/snpcalling_wild/12.Vcftools_filter_no_outgroup/indel
cd /home/vensin/workspace/snpcalling_wild/12.Vcftools_filter_no_outgroup/indel

# 2. 运行 VCFtools (INDEL)
# 新增: --remove-indv 参数
vcftools --gzvcf /home/data/10.gatk_variantfiltration/INDEL/indel_filtered.vcf.gz \
  --remove-indv O_DSMX \
  --remove-indv O_MXL \
  --remove-indv O_MZGH \
  --remove-indv O_WMMX \
  --remove-indv O_XYYG \
  --remove-indv O_ZNMX \
  --min-alleles 2 --max-alleles 2 \
  --minGQ 10 --minQ 30 --min-meanDP 6 \
  --max-missing 0.8 \
  --recode --recode-INFO-all \
  --stdout \
  2> 202_samples_indel.log \
  | bgzip -@ 8 > 202_samples_indel_filtered.recode.vcf.gz

# 3. 建立索引
tabix -p vcf 202_samples_indel_filtered.recode.vcf.gz

