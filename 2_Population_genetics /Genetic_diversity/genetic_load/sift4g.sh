## https://sift.bii.a-star.edu.sg/sift4g/SIFT4G_codes.html
## https://sift.bii.a-star.edu.sg/sift4g/SIFT4G_codes.html#SIFT%204G%20Annotator
## ***** https://chaimol.com/blog/BSA/SIFT4G-%E9%A2%84%E6%B5%8B%E6%B0%A8%E5%9F%BA%E9%85%B8%E5%8F%96%E4%BB%A3%E6%98%AF%E5%90%A6%E4%BC%9A%E5%BD%B1%E5%93%8D%E8%9B%8B%E7%99%BD%E8%B4%A8%E5%8A%9F%E8%83%BD/#gtfgffgffreadgtf9gene_biotype
# git clone --recursive https://github.com/rvaser/sift4g.git sift4g
# git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git scripts_to_build_SIFT_db
mkdir -p /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/uniref90
cd  /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/uniref90
aria2c -x 10 -s 10 -c ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
pigz -p 20 -d -k uniref90.fasta.gz

cd /home/vensin/software/SIFT_build/scripts_to_build_SIFT_db/test_files
cp homo_sapiens-test.txt guihua_config.txt


cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g
mkdir gene-annotation-src chr-src dbSNP

cd ~/software/SIFT_build/scripts_to_build_SIFT_db
## 修改make-SIFT-db-all.pl的第109行内容中，“ -d ”前面添加上-t XX.  XX 是你指定的进程数
perl /home/vensin/software/SIFT_build/scripts_to_build_SIFT_db/make-SIFT-db-all.pl \
    -config /home/vensin/workspace/snpcalling_wild/13.genetic_load/sift4g/guihua_config.txt
