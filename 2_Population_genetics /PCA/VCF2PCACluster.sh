# 1. 下载并进入目录
cd /home/vensin/software
git clone https://github.com/hewm2008/VCF2PCACluster.git
cd VCF2PCACluster
chmod 755 -R bin/*

# 2. 测试运行 (显示帮助信息则安装成功)
./bin/VCF2PCACluster -help

# 3. 添加环境变量
echo 'export PATH=$PATH:/home/vensin/software/VCF2PCACluster/bin' >> ~/.bashrc 
source ~/.bashrc

# 4. PCA
# 假设你的 pop 文件在上一级目录
# 提取 样本ID 和 Region 列，生成新的分组文件
mkdir -p /home/vensin/workspace/snpcalling_wild/12.population_genetics/pca/VCF2PCACluster
cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/pca/VCF2PCACluster

awk -F "\t" '{print $1, $3}' /home/vensin/workspace/snpcalling_wild/202samples.pop > sample.group.txt

# 检查一下格式 (应该是：SampleID Region)
head sample.group.txt

# 运行命令
VCF2PCACluster \
  -InVCF /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/208_samples_snp_filtered.LD.pruned.recode.vcf.gz \
  -InSampleGroup sample.group.txt \
  -OutPut 202_samples_AutoPCA \
  -PCnum 10 \
  -Threads 8
