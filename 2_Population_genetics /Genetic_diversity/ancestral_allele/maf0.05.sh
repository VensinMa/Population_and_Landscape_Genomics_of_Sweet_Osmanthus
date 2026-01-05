mkdir -p /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele_maf0.05
cd  /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele_maf0.05

#===========================      Workflow Steps / 流程步骤:
#===============      VCF Filtering & Data Preparation / VCF过滤与数据准备:

# ======================================================================================================================================================
mkdir -p /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs_0.05/prepare_est-sfs

# 1、Input Formatting for Est-SFS / Est-SFS 输入格式转换: 将vcf文件转换为est-sfs软件所需的输入文件格式
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs_0.05/prepare_est-sfs
wget https://raw.githubusercontent.com/VensinMa/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/refs/heads/main/2_Population_genetics%20/Genetic_diversity/ancestral_allele/vcf_to_estsfs.py

chmod +x vcf_to_estsfs.py

python vcf_to_estsfs.py  /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.maf0.05.recode.vcf  O_WMMX O_ZNMX O_MZGH
## Total sites processed: 1111139
## Successfully kept sites: 874030

# 2、Run Est-SFS / 运行 Est-SFS:
## 安装 est-sfs 
## sudo apt install libgsl-dev
## cd /home/vensin/software && wget https://github.com/VensinMa/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/raw/refs/heads/main/2_Population_genetics%20/Genetic_diversity/ancestral_allele/est-sfs-release-2.04.tar.gz
## tar -zxvf est-sfs-release-2.04.tar.gz && cd est-sfs-release-2.04
## make
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs_0.05
echo -e "n_outgroup 3\nmodel 2\nnrandom 100" > config-3outgroup.txt
echo -e "2025" > seedfile.txt
est-sfs config-3outgroup.txt  prepare_est-sfs/205_samples_snp_filtered.nomissing.maf0.05.recode_estsfs_input.txt  seedfile.txt 3_outgroup_output_file_sfs.txt  3_outgroup_output_file_p_anc.txt

# 3、VCF Polarization / VCF 极性化: 极性化原vcf文件，剔除祖先状态不确定的位点。
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs_0.05
python vcf_polarize.py /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.maf0.05.recode.vcf  \
        3_outgroup_output_file_p_anc.txt prepare_est-sfs/205_samples_snp_filtered.nomissing.maf0.05.recode_estsfs.positions.txt
## Loading data files...
## Loaded 874030 sites with ancestral probabilities.
## VCF file has been polarized based on ancestral states.
## Total sites processed in VCF: 1111139
## Sites removed (no ancestral info): 237109
## Successfully polarized sites: 874030

# 4、Genetic Load Calculation / 遗传负荷计算: 计算内类群个体突变的基因型数量 （相对于新的参考基因型 —— 祖先等位基因）
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/sift_results
python /home/vensin/software/script/indv_GT_stats_v2.py \
        /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/sift_results/205_samples_SIFT_Deleterious.vcf \
        --output  /home/vensin/workspace/snpcalling_wild/13.genetic_load/205_samples_indv_GT_stats_res_SIFT_Deleterious.txt
        
python /home/vensin/software/script/indv_GT_stats_v2.py \
        /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.snpeff_LOF.vcf \
        --output  /home/vensin/workspace/snpcalling_wild/13.genetic_load/205_samples_indv_GT_stats_res_snpeff_LOF.txt

python /home/vensin/software/script/indv_GT_stats_v2.py \
        /home/vensin/workspace/snpcalling_wild/13.genetic_load/radical/205_samples_snp_filtered.nomissing.recode_polarized.annovar_Radical.recode.vcf \
        --output  /home/vensin/workspace/snpcalling_wild/13.genetic_load/205_samples_indv_GT_stats_res_annovar_Radical.txt
