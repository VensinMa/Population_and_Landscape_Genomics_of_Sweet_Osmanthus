################################################################ SNP ################################################################
# 1. 创建输出目录
mkdir -p /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp
cd /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp

# 2. 运行 VCFtools (SNP)
# 压缩并建立索引 (生成 .vcf.gz 和 .tbi)
vcftools --gzvcf /home/data/10.gatk_variantfiltration/SNP/snp_filtered.vcf.gz \
  --remove-indv O_DSMX \
  --remove-indv O_MXL \
  --remove-indv O_MZGH \
  --remove-indv O_WMMX \
  --remove-indv O_XYYG \
  --remove-indv O_ZNMX \
  --min-alleles 2 --max-alleles 2 \
  --minGQ 10 --minQ 30 --min-meanDP 6 \
  --max-missing 0.8 --maf 0.05 \
  --recode --recode-INFO-all \
  --stdout \
  2> 202_samples_snp_filtered.recode.vcf.log \
  | bgzip -@ 8 > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz

# 然后建立索引
tabix -p vcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz



################################################################ INDEL ##############################################################
# 1. 创建并进入目录
mkdir -p /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel
cd /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel

# 2. 运行 VCFtools (INDEL) -> 管道 -> bgzip 多线程压缩
# 输入文件路径假设为: .../INDEL/indel_filtered.vcf.gz
vcftools --gzvcf /home/data/10.gatk_variantfiltration/INDEL/indel_filtered.vcf.gz \
  --remove-indv O_DSMX \
  --remove-indv O_MXL \
  --remove-indv O_MZGH \
  --remove-indv O_WMMX \
  --remove-indv O_XYYG \
  --remove-indv O_ZNMX \
  --min-alleles 2 --max-alleles 2 \
  --minGQ 10 --minQ 30 --min-meanDP 6 \
  --max-missing 0.8 --maf 0.05 \
  --recode --recode-INFO-all \
  --stdout \
  2> 202_samples_indel_filtered.recode.vcf.log \
  | bgzip -@ 8 > 202_samples_indel_filtered.recode.vcf.gz

# 3. 建立索引
tabix -p vcf 202_samples_indel_filtered.recode.vcf.gz
