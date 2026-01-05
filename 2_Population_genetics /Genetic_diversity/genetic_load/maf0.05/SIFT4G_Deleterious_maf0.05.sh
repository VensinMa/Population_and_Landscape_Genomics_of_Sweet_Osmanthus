## https://sift.bii.a-star.edu.sg/sift4g/SIFT4G_codes.html
## https://sift.bii.a-star.edu.sg/sift4g/SIFT4G_codes.html#SIFT%204G%20Annotator
## ***** https://chaimol.com/blog/BSA/SIFT4G-%E9%A2%84%E6%B5%8B%E6%B0%A8%E5%9F%BA%E9%85%B8%E5%8F%96%E4%BB%A3%E6%98%AF%E5%90%A6%E4%BC%9A%E5%BD%B1%E5%93%8D%E8%9B%8B%E7%99%BD%E8%B4%A8%E5%8A%9F%E8%83%BD/#gtfgffgffreadgtf9gene_biotype
# git clone --recursive https://github.com/rvaser/sift4g.git sift4g
# git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git scripts_to_build_SIFT_db

#=================================================== 从本地基因组和基因注释文件（.gtf）创建SIFT数据库
'''
mkdir -p /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/uniref90
cd  /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/uniref90
aria2c -x 10 -s 10 -c ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
pigz -p 20 -d -k uniref90.fasta.gz

cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g
mkdir gene-annotation-src chr-src dbSNP

cd ~/software/SIFT_build/scripts_to_build_SIFT_db
## 修改make-SIFT-db-all.pl的第109行内容中，“ -d ”前面添加上-t XX.  XX 是你指定的进程数
perl /home/vensin/software/SIFT_build/scripts_to_build_SIFT_db/make-SIFT-db-all.pl \
    -config /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/guihua_config.txt
    
## 检查数据库：数据库存储在配置文件中设置的<PARENT_DIR>/<ORG_VERSION>目录中
## The last line summarizes predictions for the entire genome: ALL # (#/#) # (#/#) # (#/#). Your database is done if the percentages are high for the last 3 different columns. Woohoo!

cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/SFZ.A && tail -n 1 CHECK_GENES.LOG
## ALL     97 (27137/27926)        99 (80301632/80758473)  86(69292747/80301632)

#=================================================== 注释 vcf
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/SFZ.A

# 2. 批量重命名 (将 Chr01.gz -> 01.gz)
# 逻辑：找到所有以 Chr 开头的文件，把 "Chr" 这三个字符删掉
for file in Chr*; do
    # 检查文件是否存在
    [ -e "$file" ] || continue
    # 执行重命名
    mv "$file" "${file#Chr}"
done

# 3. 检查结果
# 现在的列表应该是 01.gz, 02.gz, ...
ls -lh | head

cd /home/vensin/software/SIFT_build/scripts_to_build_SIFT_db  && wget https://github.com/pauline-ng/SIFT4G_Annotator/raw/master/SIFT4G_Annotator.jar
chmod +x SIFT4G_Annotator.jar
'''

cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/maf0.05
# 选择vcf和参考基因组sift4g库路径
# 设置路径变量
SIFT_JAR="/home/vensin/software/SIFT_build/scripts_to_build_SIFT_db/SIFT4G_Annotator.jar"
# 指向刚刚建好的数据库目录
SIFT_DB="/home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/SFZ.A"
# 输入文件
INPUT_VCF="/home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.maf0.05.recode.vcf"
# 输出目录
OUTPUT_DIR="/home/vensin/workspace/snpcalling_wild/13.genetic_load/maf0.05/sift4g/sift_results"

mkdir -p $OUTPUT_DIR

# 运行注释
java -Xmx10g -jar $SIFT_JAR \
    -c \
    -i $INPUT_VCF \
    -d $SIFT_DB \
    -r $OUTPUT_DIR \
    -t

'''
Started Running .......
Running in Multitranscripts mode

Chromosome      WithSIFT4GAnnotations   WithoutSIFT4GAnnotations        Progress
Chr21                   4695                    31827                   Completed : 1/23
Chr14                   6755                    40962                   Completed : 2/23
Chr19                   5098                    33866                   Completed : 3/23
Chr09                   6324                    35853                   Completed : 4/23
Chr15                   6138                    39329                   Completed : 5/23
Chr22                   5563                    30839                   Completed : 6/23
Chr13                   7981                    44350                   Completed : 7/23
Chr03                   5624                    37669                   Completed : 8/23
Chr07                   4130                    28891                   Completed : 9/23
Chr10                   6549                    43133                   Completed : 10/23
Chr08                   7363                    41387                   Completed : 11/23
Chr17                   5923                    37509                   Completed : 12/23
Chr16                   6950                    40097                   Completed : 13/23
Chr20                   5409                    35735                   Completed : 14/23
Chr23                   4995                    30613                   Completed : 15/23
Chr11                   6698                    42637                   Completed : 16/23
Chr02                   9332                    56702                   Completed : 17/23
Chr06                   8659                    45826                   Completed : 18/23
Chr05                   7703                    46941                   Completed : 19/23
Chr01                   11710                   72626                   Completed : 20/23
Chr18                   6249                    38875                   Completed : 21/23
Chr12                   7397                    43296                   Completed : 22/23
Chr04                   9090                    55841                   Completed : 23/23

Merging temp files....
SIFT4G Annotation completed !
Output directory:/home/vensin/workspace/snpcalling_wild/13.genetic_load/maf0.05/sift4g/sift_results
End Time for parallel code: Tue Jan 06 11:21:27 CST 2026
'''

#=============================== 提取 SIFT_Deleterious 位点
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/maf0.05/sift4g/sift_results
INPUT_VCF="205_samples_snp_filtered.nomissing.maf0.05.recode_SIFTpredictions.vcf"
OUTPUT_VCF="205_samples_maf0.05_SIFT_Deleterious.vcf"

# 1. 提取 VCF 表头 (Header)
grep "^#" $INPUT_VCF > $OUTPUT_VCF

# 2. 提取含有 "DELETERIOUS" 标记的行
grep "DELETERIOUS" $INPUT_VCF >> $OUTPUT_VCF

# 3. 统计一下提取到了多少个有害变异
echo "提取完成！有害变异数量："
grep -v "^#" $OUTPUT_VCF | wc -l
## 提取完成！有害变异数量：11442
