mkdir -p /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele
cd  /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele

#===========================      Workflow Steps / 流程步骤:
#===============      VCF Filtering & Data Preparation / VCF过滤与数据准备:
###  剔除亲缘关系较远的外类群，保留用于祖先状态推断的特定外类群（O_ZNMX, O_WMMX, O_MZGH）。应用严格的质量过滤标准（如无缺失数据、最小等位基因频率等）以确保位点的高质量
###  Removed distant outgroups (O_MXL, O_XYYG, O_DSMX) but retained specific outgroups (O_ZNMX, O_WMMX, O_MZGH) required for ancestral state inference. Applied strict quality filters (--max-missing 1, MA 0.01, etc.) to ensure high-quality sites without missing data.
vcftools --gzvcf /home/data/10.gatk_variantfiltration/SNP/snp_filtered.vcf.gz \
  --remove-indv O_MXL \
  --remove-indv O_XYYG \
  --remove-indv O_DSMX \
  --recode --recode-INFO-all \
  --stdout \
  2> /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.vcf.log \
  | bgzip -@ 8 > /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.vcf.gz

##  有害突变通常低频（Low frequency）稀有（Rare） (Deleterious alleles are typically rare)
##  使用 --maf 0.05（剔除频率低于 5% 的位点），会剔除那些最有可能是群体内有害突变的位点，导致严重低估群体的遗传负荷
##  --maf 0.01 是一个经验性的折衷值，足以过滤掉绝大多数随机产生的测序错误（通常是单体或双体），同时保留了那些确实存在于群体中、但受到选择压力而保持低频的真实生物学变异（有害）
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

# 1、Input Formatting for Est-SFS / Est-SFS 输入格式转换: 将vcf文件转换为est-sfs软件所需的输入文件格式
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/prepare_est-sfs
wget https://raw.githubusercontent.com/VensinMa/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/refs/heads/main/2_Population_genetics%20/Genetic_diversity/ancestral_allele/vcf_to_estsfs.py
chmod +x vcf_to_estsfs.py

python vcf_to_estsfs.py  /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.recode.vcf  O_WMMX O_ZNMX O_MZGH

# 2、Run Est-SFS / 运行 Est-SFS:
## 安装 est-sfs 
## sudo apt install libgsl-dev
## cd /home/vensin/software && wget https://github.com/VensinMa/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/raw/refs/heads/main/2_Population_genetics%20/Genetic_diversity/ancestral_allele/est-sfs-release-2.04.tar.gz
## tar -zxvf est-sfs-release-2.04.tar.gz && cd est-sfs-release-2.04
## make
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs
echo -e "n_outgroup 3\nmodel 2\nnrandom 100" > config-3outgroup.txt
echo -e "2025" > seedfile.txt
est-sfs config-3outgroup.txt  prepare_est-sfs/205_samples_snp_filtered.nomissing.recode_estsfs_input.txt  seedfile.txt 3_outgroup_output_file_sfs.txt  3_outgroup_output_file_p_anc.txt

# 3、VCF Polarization / VCF 极性化: 极性化原vcf文件，剔除祖先状态不确定的位点。
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs
python vcf_polarize.py /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.recode.vcf  \
        3_outgroup_output_file_p_anc.txt prepare_est-sfs/205_samples_snp_filtered.nomissing.recode_estsfs.positions.txt
## Loading data files...
## Loaded 3254416 sites with ancestral probabilities.
## VCF file has been polarized based on ancestral states.
## Total sites processed in VCF: 3851435
## Sites removed (no ancestral info): 597019
## Successfully polarized sites: 3254416

# 4、Genetic Load Calculation / 遗传负荷计算: 计算内类群个体突变的基因型数量 （相对于新的参考基因型 —— 祖先等位基因）
python indv_GT_stats_v2.py 205_samples_snp_filtered.nomissing.recode_polarized.vcf  --output  205_samples_snp_filtered.nomissing.recode_polarized_3_outgroup_indv_GT_stats_res.txt
## Processing 205_samples_snp_filtered.nomissing.recode_polarized.vcf...
## Writing results to 205_samples_snp_filtered.nomissing.recode_polarized_3_outgroup_indv_GT_stats_res.txt...
## Done.

