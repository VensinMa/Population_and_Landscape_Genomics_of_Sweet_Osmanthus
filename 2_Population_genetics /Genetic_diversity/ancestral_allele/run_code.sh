mkdir -p /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele
cd  /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele

###  剔除外类群样本   O_DSMX  O_MXL   O_XYYG   保留浙南木犀 网脉木犀 蒙自桂花 O_ZNMX  O_WMMX   O_MZGH
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
  --max-missing 1 --maf 0.01 \
  --recode --recode-INFO-all \
  --out /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing
## After filtering, kept 205 out of 205 Individuals
## Outputting VCF file...
## After filtering, kept 3851435 out of a possible 121081777 Sites
## Run Time = 4064.00 seconds
## rm /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.vcf.gz

# ======================================================================================================================================================
mkdir -p /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/prepare_est-sfs

# 1、将vcf文件转换为est-sfs软件所需的输入文件格式
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/prepare_est-sfs
wget https://raw.githubusercontent.com/VensinMa/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/refs/heads/main/2_Population_genetics%20/Genetic_diversity/ancestral_allele/vcf_to_estsfs.py
chmod +x vcf_to_estsfs.py

python vcf_to_estsfs.py  /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.recode.vcf  O_WMMX O_ZNMX O_MZGH

# 2、运行 est-sfs
## 安装 est-sfs 
## sudo apt install libgsl-dev
## cd /home/vensin/software && wget https://github.com/VensinMa/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/raw/refs/heads/main/2_Population_genetics%20/Genetic_diversity/ancestral_allele/est-sfs-release-2.04.tar.gz
## tar -zxvf est-sfs-release-2.04.tar.gz && cd est-sfs-release-2.04
## make
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs
echo -e "n_outgroup 3\nmodel 2\nnrandom 100" > config-3outgroup.txt
echo -e "2025" > seedfile.txt
est-sfs config-3outgroup.txt  prepare_est-sfs/205_samples_snp_filtered.nomissing.recode_estsfs_input.txt  seedfile.txt 3_outgroup_output_file_sfs.txt  3_outgroup_output_file_p_anc.txt

# 3、极性化原vcf文件
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs
python vcf_polarize.py /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.recode.vcf  \
        3_outgroup_output_file_p_anc.txt prepare_est-sfs/205_samples_snp_filtered.nomissing.recode_estsfs.positions.txt
## Loading data files...
## Loaded 3254416 sites with ancestral probabilities.
## VCF file has been polarized based on ancestral states.
## Total sites processed in VCF: 3851435
## Sites removed (no ancestral info): 597019
## Successfully polarized sites: 3254416

# 4、计算内类群个体突变的基因型数量 （相对于新的参考基因型 —— 祖先等位基因）
python indv_GT_stats_v2.py 205_samples_snp_filtered.nomissing.recode_polarized.vcf  --output  205_samples_snp_filtered.nomissing.recode_polarized_3_outgroup_indv_GT_stats_res.txt
## Processing 205_samples_snp_filtered.nomissing.recode_polarized.vcf...
## Writing results to 205_samples_snp_filtered.nomissing.recode_polarized_3_outgroup_indv_GT_stats_res.txt...
## Done.


