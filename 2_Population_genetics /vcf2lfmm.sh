##################################### ALL  ######################################
# Step 1: VCF -> .raw, VCF 转 Raw(0, 1, 2 的数字编码)
plink --vcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz \
      --recodeA --allow-extra-chr \
      --out /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.plink.recodeA

# Step 2: Raw -> LFMM
sed '1d; s/NA/9/g' /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.plink.recodeA.lfmm


##################################### LD  ######################################
# Step 1: VCF -> .raw, VCF 转 Raw(0, 1, 2 的数字编码)
plink --vcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.recode.vcf.gz \
      --recodeA --allow-extra-chr \
      --out /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.plink.recodeA

# Step 2: Raw -> LFMM
sed '1d; s/NA/9/g' /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.plink.recodeA.lfmm
