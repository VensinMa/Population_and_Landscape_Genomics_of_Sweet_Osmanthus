cd  /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele

vcftools --vcf  /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.recode.vcf \
        --max-missing 1 --maf 0.08 \
        --recode --recode-INFO-all \
        --out /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.maf0.05

## After filtering, kept 1111139 out of a possible 3851435 Sites
## Run Time = 256.00 seconds

# 提取 ID 到临时文件
grep -v "^#" /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.maf0.05.recode.vcf | awk '{print $1"\t"$2}' > target_ids.txt

# 使用 VCFtools 根据 ID 提取
vcftools --vcf /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf \
         --positions target_ids.txt \
         --recode --recode-INFO-all \
         --out /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.maf0.05.polarized
## After filtering, kept 874030 out of a possible 3254416 Sites
## Run Time = 184.00 seconds
## 205_samples_snp_filtered.nomissing.maf0.05.polarized.recode.vcf

#============================================================    
#=============================== 提取 Annovar_Radical 位点
vcftools \
    --vcf /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.maf0.05.polarized.recode.vcf \
    --positions /home/vensin/workspace/snpcalling_wild/13.genetic_load/radical/radical_positions.txt \
    --recode --recode-INFO-all \
    --out /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.maf0.05_polarized.annovar_Radical
## After filtering, kept 1997 out of a possible 874030 Sites
## Run Time = 2.00 seconds
## 2171
#============================================================
#=============================== 提取 SIFT_Deleterious 位点
vcftools --vcf /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/sift_results/205_samples_SIFT_Deleterious.vcf \
         --positions target_ids.txt \
         --recode --recode-INFO-all \
         --out /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.maf0.05.polarized.SIFT_Deleterious
## After filtering, kept 8777 out of a possible 56051 Sites
## Run Time = 2.00 seconds
## 14031

#============================================================
#=============================== 提取 SNPEFF_LOF 位点
vcftools --vcf /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.snpeff_LOF.vcf \
         --positions target_ids.txt \
         --recode --recode-INFO-all \
         --out /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.maf0.05.polarized.SNPEFF_LOF
## After filtering, kept 514 out of a possible 3050 Sites
## Run Time = 0.00 seconds
## 825


#============================================================
python /home/vensin/software/script/indv_GT_stats_v2.py \
        205_samples_snp_filtered.nomissing.maf0.05_polarized.annovar_Radical.recode.vcf \
        --output  205_sample_maf0.05s_indv_GT_stats_res_SIFT_Deleterious.txt
        
python /home/vensin/software/script/indv_GT_stats_v2.py \
        205_samples_snp_filtered.nomissing.maf0.05.polarized.SIFT_Deleterious.recode.vcf \
        --output  205_samples_maf0.05_indv_GT_stats_res_snpeff_LOF.txt

python /home/vensin/software/script/indv_GT_stats_v2.py \
        205_samples_snp_filtered.nomissing.maf0.05.polarized.SNPEFF_LOF.recode.vcf \
        --output  205_samples_maf0.05_indv_GT_stats_res_annovar_Radical.txt
