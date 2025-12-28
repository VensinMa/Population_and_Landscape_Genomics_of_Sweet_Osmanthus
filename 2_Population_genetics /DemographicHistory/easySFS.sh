cd /home/vensin/software
conda create -n easySFS && conda activate easySFS
git clone https://github.com/isaacovercast/easySFS.git

cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/easySFS
easySFS.py -i /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf -p 205samples.pop -a --unfolded --preview
# -a 表示生成 Unfolded SFS (适用于已极化的数据)
# --proj 顺序必须严格对应: East, Central, SW-Yunnan, SW-Guizhou, OUTGROUP

easySFS.py -i /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf \
-p 205samples.pop \
-a \
--unfolded\
--proj 272,54,24,54,6 \
-o output_fitcoal_unfolded

######################################################################## 
mkdir -p  /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/easySFS/135samples
cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/easySFS/135samples
easySFS.py -i /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf -p ../135samples.pop -a --unfolded --preview
# -a 表示生成 Unfolded SFS (适用于已极化的数据)
'''
(base) vensin@Ubuntu24:~/workspace/snpcalling_wild/12.population_genetics/Demographic_History/easySFS/135samples$ easySFS.py -i /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf -p ../135samples.pop -a --unfolded --preview
  Processing 4 populations - ['East', 'Central', 'Southwest-Yunnan', 'Southwest-Guizhou']
Samples in pops file not present in VCF: 
Samples in VCF not present in pops file: HJLL_2, HJLL_3, O_MZGH, LX_4, JMX_16, XC_2, ZJP_3, HYX_2, O_ZNMX, ZJP_10, YX_7, FZA_6, XC_5, YZY_10, YZY_4, QDH_1, HYX_6, HJLL_7, GHX_13, FZA_1, HYX_3, HJLL_6, DST_9, ZJP_1, ST_1, XC_1, ZJS_10, HJLL_1, GHX_4, GHX_1, ZJP_7, EJ_6, ST_6, FZA_4, QDH_3, JMX_3, SL_2, JMX_2, ST_5, ZJP_8, EJ_7, JMX_15, GHX_2, ZJP_5, LX_2, ST_4, ST_3, ZJP_2, GHX_10, XC_6, O_WMMX, HYX_4, DST_5, EJ_1, FZA_5, EJ_8, FZA_3, YZY_2, ZJP_4, GHX_9, GHX_8, FZA_8, SL_3, ST_2, HYX_1, EJ_3, ST_8, ST_7, HYX_5, EJ_4
Continue, excluding samples not in both pops file and VCF? (yes/no)
yes
'''

# --proj 顺序必须严格对应: East, Central, SW-Yunnan, SW-Guizhou, OUTGROUP
easySFS.py -i /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf \
-p ../135samples.pop \
-a \
--unfolded \
--proj 202,22,24,22 \
-o output_fitcoal_unfolded_135samples


'''
Processing 4 populations - ['East', 'Central', 'Southwest-Yunnan', 'Southwest-Guizhou']
Samples in pops file not present in VCF: 
Samples in VCF not present in pops file: ZJP_1, HYX_2, EJ_7, DST_5, GHX_4, ST_6, YZY_2, XC_1, ST_2, EJ_3, ZJP_3, HJLL_6, ZJP_5, GHX_1, FZA_3, HYX_3, YZY_4, JMX_16, JMX_2, ZJP_2, ST_5, QDH_1, EJ_8, HJLL_3, QDH_3, HJLL_2, FZA_1, ST_4, LX_2, O_ZNMX, ST_8, HYX_6, GHX_2, JMX_15, DST_9, EJ_6, YX_7, O_WMMX, GHX_13, JMX_3, ZJP_8, ZJP_4, HJLL_7, ST_1, EJ_1, FZA_5, GHX_8, FZA_8, ZJP_10, GHX_9, ST_7, HJLL_1, XC_5, HYX_4, FZA_4, GHX_10, YZY_10, FZA_6, ZJP_7, XC_6, SL_2, O_MZGH, ST_3, HYX_1, SL_3, EJ_4, HYX_5, XC_2, LX_4, ZJS_10
Continue, excluding samples not in both pops file and VCF? (yes/no)
yes
SFS files written to /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/easySFS/135samples/output_fitcoal_unfolded_135samples
'''
