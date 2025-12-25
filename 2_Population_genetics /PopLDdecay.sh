conda create -n poplddecay
conda activate poplddecay
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install -c bioconda -c conda-forge poplddecay
echo 'export PATH="~/anaconda3/envs/poplddecay/bin:$PATH"' >> ~/.bashrc



cd  /home/vensin/workspace/snpcalling_wild/12.population_genetics/PopLDdecay

### Step1: 计算 LD decay
PopLDdecay  -InVCF /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz -MaxDist 300 -SubPop East.txt -OutStat East.PopLDdecay
PopLDdecay  -InVCF /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz -MaxDist 300 -SubPop Central.txt -OutStat Central.PopLDdecay
PopLDdecay  -InVCF /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz -MaxDist 300 -SubPop Southwest-Guizhou.txt -OutStat Southwest-GZ.PopLDdecay
PopLDdecay  -InVCF /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.recode.vcf.gz -MaxDist 300 -SubPop Southwest-Yunnan.txt -OutStat Southwest-YN.PopLDdecay


### Step2: 绘图
## 单个谱系绘图
Plot_OnePop.pl -inFile East.PopLDdecay.stat.gz --output East.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile Central.PopLDdecay.stat.gz --output Central.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile Southwest-GZ.PopLDdecay.stat.gz --output Southwest-GZ.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile Southwest-YN.PopLDdecay.stat.gz --output Southwest-YN.PopLDdecay.stat -bin1 10 -bin2 100



## 多谱系绘图
cat <<EOF > K4_multi.list
East.PopLDdecay.stat.gz East
Central.PopLDdecay.stat.gz Central
Southwest-GZ.PopLDdecay.stat.gz Southwest-GZ
Southwest-YN.PopLDdecay.stat.gz Southwest-YN
EOF


Plot_MultiPop.pl -inList K4_multi.list --output K4.PopLDdecay -bin1 10 -bin2 100
Plot_MultiPop.pl -inList K2_multi.list --output K2.PopLDdecay -bin1 10 -bin2 100

gzip -d -k Central-East.PopLDdecay.stat.bin.gz Central.PopLDdecay.stat.bin.gz  East.PopLDdecay.stat.bin.gz  West-GZ.PopLDdecay.stat.bin.gz  West-YN.PopLDdecay.stat.bin.gz  West.PopLDdecay.stat.bin.gz
