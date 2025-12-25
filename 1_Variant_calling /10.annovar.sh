##  annovar官方网站的的下载需要使用邮箱注册后才可下载
# https://www.openbioinformatics.org/annovar/annovar_download_form.php

################################# annovar 安装 ##############################################
cd /home/vensin/software/ && wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz  
tar -zxvf annovar.latest.tar.gz

## 安装注释文件格式转换工具 gffread gff3ToGenePred gtfToGenePred
conda create -n ucsc
conda activate ucsc
conda install bioconda::ucsc-gff3togenepred
conda install bioconda::ucsc-gtftogenepred

## 添加环境变量 export PATH
echo 'export PATH=/home/vensin/software/annovar:$PATH' >> ~/.bashrc
echo 'export PATH=/home/vensin/anaconda3/envs/ucsc/bin/:$PATH' >> ~/.bashrc
source ~/.bashrc

########################  gff格式转gtf格式  ############################################################################
mkdir -p /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/Annovar/result && cd /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/Annovar/
gffread  /home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.gff -T -o SFZ.A.onlychr.gtf
gtfToGenePred  -genePredExt SFZ.A.onlychr.gtf SFZ.A.onlychr_refGene.txt

retrieve_seq_from_fasta.pl --format refGene --seqfile /home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa  SFZ.A.onlychr_refGene.txt --out SFZ.A.onlychr_refGeneMrna.fa

#############################  生成表格格式输入文件  ###################################################################
convert2annovar.pl -format vcf4 -allsample -withfreq \
  /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz \
  > 202_samples_snp_filtered.recode.annovar.input

## 进行变异注释 (如果需要所有信息，添加-separate) -separate 将每种类型的变异分开注释到不同的文件中
annotate_variation.pl -geneanno --neargene 2000 -buildver  SFZ.A.onlychr -dbtype refGene -outfile SFZ.A.onlychr.snp.annovar -exonsort 202_samples_snp_filtered.recode.annovar.input  ./

####################################   统计每种类型（突变位置）SNP的数量  ###############################################
# cat LYG.hic.snp.annovar.variant_function | cut -f 1 | sed 's/;/\n/g' | sort | uniq -c
cat SFZ.A.onlychr.snp.annovar.variant_function | cut -f 1 | sort | uniq -c

1434710 downstream
 710456 exonic
      5 exonic;splicing
8028563 intergenic
2929085 intronic
   3718 splicing
1731234 upstream
 269330 upstream;downstream
 232506 UTR3
 152326 UTR5
   1741 UTR5;UTR3
#######################################   统计外显子区域不同突变类型SNP的数量  #############################################
cat SFZ.A.onlychr.snp.annovar.exonic_variant_function | awk '{print $2}' | sort | uniq -c

390954 nonsynonymous
   8859 stopgain
   1132 stoploss
 309516 synonymous

 
