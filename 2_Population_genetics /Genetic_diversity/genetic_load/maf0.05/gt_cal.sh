#============================================================
python /home/vensin/software/script/indv_GT_stats_v2.py \
         /home/vensin/workspace/snpcalling_wild/13.genetic_load/maf0.05/205_samples_snp_filtered.nomissing.maf0.05.recode_polarized.annovar_Radical.recode.vcf \
        --output  205_sample_maf0.05s_indv_GT_stats_res_SIFT_Deleterious.txt
        
python /home/vensin/software/script/indv_GT_stats_v2.py \
        /home/vensin/workspace/snpcalling_wild/13.genetic_load/maf0.05/sift4g/sift_results/205_samples_maf0.05_SIFT_Deleterious.vcf \
        --output  205_samples_maf0.05_indv_GT_stats_res_snpeff_LOF.txt

python /home/vensin/software/script/indv_GT_stats_v2.py \
        /home/vensin/workspace/snpcalling_wild/13.genetic_load/maf0.05/205_samples_snp_filtered.nomissing.maf0.05_polarized.snpeff_LOF.vcf \
        --output  205_samples_maf0.05_indv_GT_stats_res_annovar_Radical.txt
