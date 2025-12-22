mkdir -p /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele
cd  /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele

###  剔除外类群样本   O_DSMX  O_MXL   O_XYWJM   保留浙南木犀 网脉木犀 蒙自桂花 O_ZNMX  NFM_1   O_MZGH
vcftools --gzvcf /home/data/10.gatk_variantfiltration/SNP/snp_filtered.vcf.gz \
  --remove-indv O_MXL \
  --remove-indv O_XYYG \
  --remove-indv O_DSMX \
  --recode --recode-INFO-all \
  --stdout \
  2> /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.vcf.log \
  | bgzip -@ 8 > /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.vcf.gz

vcftools --gzvcf /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.vcf.gz \
  --min-alleles 2 --max-alleles 2 \
  --minGQ 10 --minQ 30 --min-meanDP 6 \
  --max-missing 1 --maf 0.05 \
  --recode --recode-INFO-all \
  --out /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing
rm /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.vcf.gz
# ======================================================================================================================================================
mkdir -p /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/prepare_est-sfs
cd  /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/prepare_est-sfs

plink --vcf /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.recode.vcf \
      --out 205_samples_snp_filtered.nomissing.recode.plink \
      --recode vcf-iid  --allow-extra-chr  --keep-allele-order  --set-missing-var-ids @:# 


# 1、将vcf文件转换为est-sfs软件所需的输入文件格式
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/prepare_est-sfs
python vcf_to_estsfs.py  205_samples_snp_filtered.nomissing.recode.plink  O_WMMX O_ZNMX O_MZGH

'''
(base) vensin@ubuntu24-04:~/workspace/est-sfs/prepare_est-sfs$ python vcf_to_estsfs.py  196samples_filtered_2_outgroup.snp.nomissing.rename.plink.vcf  O-WMMX O-ZNMX
VCF file has been converted to est-sfs input file.
Total sites processed: 1104492
Successfully kept sites: 948992
(base) vensin@ubuntu24-04:~/workspace/est-sfs/prepare_est-sfs$ python vcf_to_estsfs.py  197samples_filtered_3_outgroup.snp.nomissing.rename.plink.vcf  O-WMMX O-ZNMX O-MZGH
VCF file has been converted to est-sfs input file.
Total sites processed: 996932
Successfully kept sites: 808574
'''
'''
-rw-rw-r-- 1 vensin vensin  15M 11月 18 19:20 196samples_filtered_2_outgroup.snp.nomissing.rename.plink_processing.log
-rw-rw-r-- 1 vensin vensin  22M 11月 18 19:20 196samples_filtered_2_outgroup.snp.nomissing.rename.plink_estsfs.positions.txt
-rw-rw-r-- 1 vensin vensin  25M 11月 18 19:20 196samples_filtered_2_outgroup.snp.nomissing.rename.plink_estsfs_input.txt
drwxrwxr-x 2 vensin vensin 4.0K 11月 18 19:20 ./
-rw-rw-r-- 1 vensin vensin  18M 11月 18 19:20 197samples_filtered_3_outgroup.snp.nomissing.rename.plink_processing.log
-rw-rw-r-- 1 vensin vensin  19M 11月 18 19:20 197samples_filtered_3_outgroup.snp.nomissing.rename.plink_estsfs.positions.txt
-rw-rw-r-- 1 vensin vensin  28M 11月 18 19:20 197samples_filtered_3_outgroup.snp.nomissing.rename.plink_estsfs_input.txt
'''

# 2、运行 est-sfs
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/
'''
(base) vensin@ubuntu24-04:~/workspace/est-sfs$ cat config-2outgroup.txt config-3outgroup.txt
n_outgroup 2
model 2
nrandom 100
n_outgroup 3
model 2
nrandom 100
'''
est-sfs config-3outgroup.txt  /home/vensin/workspace/est-sfs/prepare_est-sfs/197samples_filtered_3_outgroup.snp.nomissing.rename.plink_estsfs_input.txt  seedfile.txt 3_outgroup_output_file_sfs.txt  3_outgroup_output_file_p_anc.txt

# 3、极性化原vcf文件
cd /home/vensin/workspace/est-sfs/
python vcf_polarize.py /home/vensin/workspace/est-sfs/prepare_est-sfs/196samples_filtered_2_outgroup.snp.nomissing.rename.plink.vcf  2_outgroup_output_file_p_anc.txt /home/vensin/workspace/est-sfs/prepare_est-sfs/196samples_filtered_2_outgroup.snp.nomissing.rename.plink_estsfs.positions.txt
# Total polarized sites: 948992
python vcf_polarize.py /home/vensin/workspace/est-sfs/prepare_est-sfs/197samples_filtered_3_outgroup.snp.nomissing.rename.plink.vcf  3_outgroup_output_file_p_anc.txt /home/vensin/workspace/est-sfs/prepare_est-sfs/197samples_filtered_3_outgroup.snp.nomissing.rename.plink_estsfs.positions.txt
# Total polarized sites: 808574

# 4、计算内类群个体突变的基因型数量 （相对于新的参考基因型 —— 祖先等位基因）
python indv_GT_stats.py 196samples_filtered_2_outgroup.snp.nomissing.rename.plink.polarized.vcf  --output  196samples_filtered_2_outgroup_indv_GT_stats_res.txt
python indv_GT_stats.py 197samples_filtered_3_outgroup.snp.nomissing.rename.plink.polarized.vcf  --output  197samples_filtered_3_outgroup_indv_GT_stats_res.txt
