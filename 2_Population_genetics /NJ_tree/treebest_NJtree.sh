cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/tree_treebest
python /home/vensin/software/vcf2phylip-2.8/vcf2phylip.py \
      --input /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/208_samples_snp_filtered.LD.pruned.recode.vcf.gz \
      --fasta --output-prefix 208_samples_snp_filtered.LD.pruned.treebest
      
treebest nj -b 1000 208_samples_snp_filtered.LD.pruned.treebest.min4.fasta > 208_samples_snp_filtered.LD.pruned.treebest.out
