##################################### ALL SNP  ######################################
# Step 1: VCF -> .raw, VCF 转 Raw(0, 1, 2 的数字编码)
plink --vcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz \
      --recodeA --allow-extra-chr \
      --out /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.plink.recodeA

# Step 2: Raw -> LFMM
sed '1d; s/NA/9/g' /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.plink.recodeA.lfmm

# 
bcftools query -f '%CHROM:%POS\n' /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.ID

##################################### ALL SNP + --max-missing 0.95 (缺失0.05)  ######################################
vcftools --gzvcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz \
  --max-missing 0.95 --maf 0.05 \
  --recode --recode-INFO-all \
  --stdout \
  2> 202_samples_snp_filtered.missing0.05.recode.vcf.log \
  | bgzip -@ 8 > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.missing0.05.recode.vcf.gz

# Step 1: VCF -> .raw, VCF 转 Raw(0, 1, 2 的数字编码)
plink --vcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.missing0.05.recode.vcf.gz \
      --recodeA --allow-extra-chr \
      --out /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.missing0.05.plink.recodeA

# Step 2: Raw -> LFMM
sed '1d; s/NA/9/g' /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.missing0.05.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.missing0.05.plink.recodeA.lfmm

# Step 3: 提取变异位点ID
bcftools query -f '%CHROM:%POS\n' /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.missing0.05.recode.vcf.gz > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.missing0.05.ID  
  
##################################### LD SNP ######################################
# Step 1: VCF -> .raw, VCF 转 Raw(0, 1, 2 的数字编码)
plink --vcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.recode.vcf.gz \
      --recodeA --allow-extra-chr \
      --out /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.plink.recodeA

# Step 2: Raw -> LFMM
sed '1d; s/NA/9/g' /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.plink.recodeA.lfmm


##################################### ALL INDEL ######################################
# Step 1: VCF -> .raw, VCF 转 Raw(0, 1, 2 的数字编码)
plink --vcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.recode.vcf.gz \
      --recodeA --allow-extra-chr \
      --out /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.plink.recodeA

# Step 2: Raw -> LFMM
sed '1d; s/NA/9/g' /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.plink.recodeA.lfmm

# Step 3: 提取变异位点ID
bcftools query -f '%CHROM:%POS\n' /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.recode.vcf.gz > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.ID

##################################### ALL INDEL + --max-missing 0.95 (缺失0.05)  ######################################
vcftools --gzvcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.recode.vcf.gz \
  --max-missing 0.95 --maf 0.05 \
  --recode --recode-INFO-all \
  --stdout \
  2> 202_samples_indel_filtered.missing0.05.recode.vcf.log \
  | bgzip -@ 8 > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.missing0.05.recode.vcf.gz
## After filtering, kept 202 out of 202 Individuals
## Outputting VCF file...
## After filtering, kept 810159 out of a possible 1163368 Sites
## Run Time = 186.00 seconds

# Step 1: VCF -> .raw, VCF 转 Raw(0, 1, 2 的数字编码)
plink --vcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.missing0.05.recode.vcf.gz \
      --recodeA --allow-extra-chr \
      --out /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.missing0.05.plink.recodeA

# Step 2: Raw -> LFMM
sed '1d; s/NA/9/g' /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.missing0.05.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.missing0.05.plink.recodeA.lfmm

# Step 3: 提取变异位点ID
bcftools query -f '%CHROM:%POS\n' /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.missing0.05.recode.vcf.gz > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.missing0.05.ID  
  



      
