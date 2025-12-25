cd /home/vensin/software
conda create -n easySFS && conda activate easySFS
git clone https://github.com/isaacovercast/easySFS.git

./easySFS.py -i /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.recode.vcf -p pops_file.txt --preview
