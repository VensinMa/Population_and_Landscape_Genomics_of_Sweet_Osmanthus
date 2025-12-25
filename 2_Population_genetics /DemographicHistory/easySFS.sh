cd /home/vensin/software
conda create -n easySFS && conda activate easySFS
git clone https://github.com/isaacovercast/easySFS.git

cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/easySFS
easySFS.py -i /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf -p 205samples.pop -a --preview
# -a 表示生成 Unfolded SFS (适用于已极化的数据)
# --proj 顺序必须严格对应: East, Central, SW-Yunnan, SW-Guizhou, OUTGROUP

easySFS.py -i /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf \
-p 205samples.pop \
-a \
--proj 272,54,24,54,6 \
-o output_fitcoal_unfolded
