## https://sift.bii.a-star.edu.sg/sift4g/SIFT4G_codes.html
## https://sift.bii.a-star.edu.sg/sift4g/SIFT4G_codes.html#SIFT%204G%20Annotator
## ***** https://chaimol.com/blog/BSA/SIFT4G-%E9%A2%84%E6%B5%8B%E6%B0%A8%E5%9F%BA%E9%85%B8%E5%8F%96%E4%BB%A3%E6%98%AF%E5%90%A6%E4%BC%9A%E5%BD%B1%E5%93%8D%E8%9B%8B%E7%99%BD%E8%B4%A8%E5%8A%9F%E8%83%BD/#gtfgffgffreadgtf9gene_biotype
# git clone --recursive https://github.com/rvaser/sift4g.git sift4g
# git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git scripts_to_build_SIFT_db

#=================================================== 从本地基因组和基因注释文件（.gtf）创建SIFT数据库

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
    
'''
entered mkdir /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/SFZ.A
converting gene format to use-able input
done converting gene format
making single records file
done making single records template
making noncoding records file
done making noncoding records
make the fasta sequences
done making the fasta sequences
start siftsharp, getting the alignments
/home/vensin/software/SIFT_build/sift4g/bin/sift4g -t 24 -d /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/uniref90/uniref90.fasta -q /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/all_prot.fasta --subst /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/subst --out /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/SIFT_predictions --sub-results 
** Checking query data and substitutions files **
* processing queries: 100.00/100.00% *

** Searching database for candidate sequences **
* processing database part 352 (size ~0.25 GB): 100.00/100.00% *

** Aligning queries with candidate sequences **
* processing database part 88 (size ~1.00 GB): 100.00/100.00% *

** Selecting alignments with median threshold: 2.75 **
* processing queries: 100.00/100.00% *

** Generating SIFT predictions with sequence identity: 100.00% **
* processing queries: 100.00/100.00% *

done getting all the scores
populating databases
checking the databases
zipping up /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/chr-src/*
All done!
'''

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

cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g
# 选择vcf和参考基因组sift4g库路径
# 设置路径变量
SIFT_JAR="/home/vensin/software/SIFT_build/scripts_to_build_SIFT_db/SIFT4G_Annotator.jar"
# 指向刚刚建好的数据库目录
SIFT_DB="/home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/SFZ.A"
# 输入文件
INPUT_VCF="/home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf"
# 输出目录
OUTPUT_DIR="/home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/sift_results"

mkdir -p $OUTPUT_DIR

# 运行注释
java -Xmx10g -jar $SIFT_JAR \
    -c \
    -i $INPUT_VCF \
    -d $SIFT_DB \
    -r $OUTPUT_DIR \
    -t
'''
Start Time for SIFT4G code: Thu Dec 25 21:02:34 CST 2025

Started Running .......
Running in Multitranscripts mode

Chromosome      WithSIFT4GAnnotations   WithoutSIFT4GAnnotations        Progress
Chr19                   13497                   98806                   Completed : 1/23
Chr21                   11888                   87502                   Completed : 2/23
Chr07                   11156                   86996                   Completed : 3/23
Chr20                   14548                   107489                  Completed : 4/23
Chr03                   15935                   109744                  Completed : 5/23
Chr23                   13972                   96144                   Completed : 6/23
Chr09                   16694                   107335                  Completed : 7/23
Chr18                   16259                   111093                  Completed : 8/23
Chr17                   16584                   117352                  Completed : 9/23
Chr22                   14110                   87037                   Completed : 10/23
Chr16                   18737                   116893                  Completed : 11/23
Chr12                   20292                   130873                  Completed : 12/23
Chr15                   18124                   122026                  Completed : 13/23
Chr08                   20396                   125396                  Completed : 14/23
Chr14                   17453                   116333                  Completed : 15/23
Chr11                   17993                   121893                  Completed : 16/23
Chr10                   17236                   126298                  Completed : 17/23
Chr13                   21725                   134781                  Completed : 18/23
Chr02                   24665                   166466                  Completed : 19/23
Chr01                   32034                   223536                  Completed : 20/23
Chr04                   22944                   159701                  Completed : 21/23
Chr05                   20608                   141075                  Completed : 22/23
Chr06                   23116                   139681                  Completed : 23/23

Merging temp files....
SIFT4G Annotation completed !
Output directory:/home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/sift_results
End Time for parallel code: Thu Dec 25 21:06:24 CST 2025
'''

#=============================== 提取 SIFT_Deleterious 位点
INPUT_VCF="205_samples_snp_filtered.nomissing.recode_polarized_SIFTpredictions.vcf"
OUTPUT_VCF="205_samples_SIFT_Deleterious.vcf"

# 1. 提取 VCF 表头 (Header)
grep "^#" $INPUT_VCF > $OUTPUT_VCF

# 2. 提取含有 "DELETERIOUS" 标记的行
grep "DELETERIOUS" $INPUT_VCF >> $OUTPUT_VCF

# 3. 统计一下提取到了多少个有害变异
echo "提取完成！有害变异数量："
grep -v "^#" $OUTPUT_VCF | wc -l
## 提取完成！有害变异数量：56051








