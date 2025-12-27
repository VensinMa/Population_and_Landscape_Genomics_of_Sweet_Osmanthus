## NR注释记录 20240325
### 1、安装 aspera 与下载 NR 蛋白数据库文件

# 下载安装 aspera 软件
# conda create -n aspera
# conda activate aspera
# conda install -c hcc aspera-cli -y

# 添加并立即生效环境变量
# echo 'export PATH="/home/vensin/anaconda3/envs/aspera/bin/:$PATH"' >> ~/.bashrc
# source ~/.bashrc

# 检查是否安装成功
# ascp -h

# 切换到工作目录 workspace 
cd /home/data/11.vcftools_filter/snp/nr.annotation
    
# 下载 nr 蛋白库
#ascp -i /home/vensin/anaconda3/envs/aspera/etc/asperaweb_id_dsa.openssh -l 1000M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/FASTA/nr.gz ./ &
#ascp -i /home/vensin/anaconda3/envs/aspera/etc/asperaweb_id_dsa.openssh -l 1000M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/FASTA/nr.gz.md5 ./ &
# wget https://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz.md5
md5sum -c nr.gz.md5
# gunzip -c nr.gz > nr.fasta
unpigz -c -p 4 nr.gz > nr.fasta


# 下载 Nr 数据库中蛋白编号 mprot.accession 和物种编号 taxid 的对应关系信息 
# ascp -i /home/vensin/anaconda3/envs/aspera/etc/asperaweb_id_dsa.openssh -l 1000M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/prot.accession2taxid.gz ./ &
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

# 下载 Nr 数据库中物种编号 taxdmp 的层次信息 
# ascp -i /home/vensin/anaconda3/envs/aspera/etc/asperaweb_id_dsa.openssh -l 1000M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip

gffread  /home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.gff  -g  /home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa  -w  SFZ.A.gene.fasta
gffread  /home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.gff  -g  /home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa  -x  SFZ.A.cds.fasta
gffread  /home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.gff  -g  /home/vensin/workspace/snpcalling_wild/0.genome/SFZ.A.onlychr.fa  -y  SFZ.A.pep.fasta

### 2、使用 Blast/Diamond 进行 NR 注释
cd /home/data/11.vcftools_filter/snp/nr.annotation
# 使用diamond  软件的子命令 makedb 将 fasta 格式的蛋白序列创建成后缀为 dmnd 的数据库文件：
diamond makedb --in nr.fasta --db nr.db --threads 28

# 将物种全基因组核酸/蛋白序列 blastx / blastp 到构建好的数据库：
diamond blastp --db nr.db.dmnd --query SFZ.A.pep.fasta \
    --out SFZ.A.pep.Nr.annotations \
    --outfmt 6 qseqid sseqid pident evalue bitscore qlen slen length mismatch gapopen qstart qend sstart send stitle \
    --sensitive --max-target-seqs 5 --evalue 1e-5 --threads 28

## Total time = 24324.1s
## Reported 213684 pairwise alignments, 213684 HSPs.
## 43374 queries aligned.

python unique_gene_annotations.py SFZ.A.pep.Nr.annotations
## 正在读取: SFZ.A.pep.Nr.annotations
## 正在写入: SFZ.A.pep.uniq.Nr.annotations ...
## 处理完成！共保留了 43374 个基因的最佳注释。
 python merge_annovar_nr_v2.py \
     /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/Annovar/SFZ.A.onlychr.snp.annovar.exonic_variant_function \
     /home/data/11.vcftools_filter/snp/nr.annotation/SFZ.A.pep.uniq.Nr.annotations \
     SFZ.A.final.snp.annovar_nr.annotation.xls
## 正在加载注释库: /home/data/11.vcftools_filter/snp/nr.annotation/SFZ.A.pep.uniq.Nr.annotations ...
## 正在处理并合并: /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/Annovar/SFZ.A.onlychr.snp.annovar.exonic_variant_function ...
## 合并完成！
## 输出文件: SFZ.A.final.snp.annovar_nr.annotation.xls
## 总处理行数: 710461
## 成功匹配注释: 690798
