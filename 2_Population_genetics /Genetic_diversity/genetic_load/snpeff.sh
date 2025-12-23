####  https://blog.csdn.net/2302_79242191/article/details/134630776?utm_medium=distribute.pc_relevant.none-task-blog-2~default~baidujs_baidulandingword~default-0-134630776-blog-105530364.235^v43^pc_blog_bottom_relevance_base2&spm=1001.2101.3001.4242.1&utm_relevant_index=1
## https://pcingola.github.io/SnpEff/snpeff/introduction/
### 0.准备
##cd /home/vensin/software/ && wget https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip
##unzip snpEff_latest_core.zip && cd /home/vensin/software/snpEff  && chmod +x SnpSift.jar  && chmod +x snpEff.jar
##cd /home/vensin/software/snpEff_latest_core/snpEff  && chmod +x SnpSift.jar  && chmod +x snpEff.jar
echo "guihua.genome : Osmanthus fragrans" >> /home/vensin/software/snpEff/snpEff.config

##  mkdir -p data/XXX  
mkdir -p /home/vensin/software/snpEff/data/guihua 

## 在data下面存放两个文件  SFZ.A.onlychr.fa: 参考基因组  SFZ.A.onlychr.gff: 注释文件，GFF3格式
cd /home/vensin/software/snpEff/data/
mv SFZ.A.onlychr.fa sequences.fa
mv SFZ.A.onlychr.gff genes.gff
## java -jar snpEff.jar build -gff3 -v XXX

### 1.建库
cd /home/vensin/software/snpEff/
java -jar /home/vensin/software/snpEff/snpEff.jar build -gff3 -v guihua -noCheckProtein -noCheckCds

### 2.snpeff 注释  
# 输入文件路径
INPUT_VCF="/home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf"
# 输出文件路径
OUTPUT_VCF="/home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.snpeff.vcf"

cd /home/vensin/software/snpEff/
java -Xmx20g -jar snpEff.jar \
    -c snpEff.config \
    -v guihua \
    -ud 2000 \
    -csvStats guihua.snpeff.205.csv \
    -htmlStats guihua.snpeff.205.html \
    $INPUT_VCF > $OUTPUT_VCF

### 3.提取注释结果中的 LOF 位点
# 使用 SnpSift filter
# 逻辑：保留 Header，并且保留 INFO 字段中存在 "LOF" 标记的行
DIR="/home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs"
INPUT="$DIR/205_samples_snp_filtered.nomissing.recode_polarized.snpeff.vcf"
OUTPUT="$DIR/205_samples_snp_filtered.nomissing.recode_polarized.snpeff_LOF.vcf"
java -jar /home/vensin/software/snpEff/SnpSift.jar filter "( exists LOF )" $INPUT > $OUTPUT
grep -v "^#" /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.snpeff_LOF.vcf | wc -l
## 3050

