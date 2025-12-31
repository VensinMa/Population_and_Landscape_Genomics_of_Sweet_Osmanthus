# 1. 进入工作目录
cd /home/vensin/workspace/snpcalling_wild/13.genetic_load/radical

# 2. 复制矩阵文件 (如果还没复制)
cp /home/vensin/software/annovar/example/grantham.matrix ./

# 3. 格式转换 (VCF -> Annovar)
# 注意：如果输入的是标准VCF，用 vcf4 可能比 vcf4old 更稳妥，但如果确定是旧格式则保持不变
convert2annovar.pl \
    --format vcf4old \
    /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf \
    --outfile 205_samples_snp_filtered.nomissing.recode_polarized.avinput

'''
(base) vensin@Ubuntu24:~/workspace/snpcalling_wild/13.genetic_load/radical$ convert2annovar.pl \
    --format vcf4old \
    /home/vensin/workspace/snpcalling_wild/13.genetic_load/est-sfs/205_samples_snp_filtered.nomissing.recode_polarized.vcf \
    --outfile 205_samples_snp_filtered.nomissing.recode_polarized.avinput
NOTICE: for SNPs, column 6 and beyond MAY BE heterozygosity status, quality score, read depth, RMS mapping quality, quality by depth, if these information can be recognized automatically
NOTICE: for indels, column 6 and beyond MAY BE heterozygosity status, quality score, read depth, read count supporting indel call, RMS mapping quality, if these information can be recognized automatically
NOTICE: Read 3254483 lines and wrote 3055651 different variants at 3254416 genomic positions (3254416 SNPs and 0 indels)
NOTICE: Among 3254416 different variants at 3254416 positions, 231007 are heterozygotes, 2824644 are homozygotes
NOTICE: Among 3254416 SNPs, 2054336 are transitions, 1200080 are transversions (ratio=1.71)
'''

# 4. 运行注释 (核心步骤)
# 数据库目录指向你存放 SFZ.A.onlychr_refGene.txt 的地方
annotate_variation.pl \
    -buildver SFZ.A.onlychr \
    -outfile SFZ.A.grantham \
    --aamatrixfile grantham.matrix \
    205_samples_snp_filtered.nomissing.recode_polarized.avinput \
    /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/Annovar/

# 5. 提取所有非同义突变
grep 'nonsynonymous' SFZ.A.grantham.exonic_variant_function > nonsynonymous_variants.txt

# 6. 筛选 Grantham Score > 150 的激进突变
# 解释：awk 根据 "AAMatrix=" 切分，取后面的数字判断是否 > 150
awk -F'AAMatrix=' '{split($2, a, ","); if(a[1] > 150) print $0}' nonsynonymous_variants.txt > SFZ.A.radical_sites_list.txt

# 7. 统计并展示结果数量
echo "--------------------------------"
echo "非同义突变总数:"
wc -l nonsynonymous_variants.txt
echo "Radical (激进) 突变数 (>150):"
wc -l SFZ.A.radical_sites_list.txt
echo "--------------------------------"

## 非同义突变总数:
## 188775 nonsynonymous_variants.txt
## Radical (激进) 突变数 (>150):
## 11702 SFZ.A.radical_sites_list.txt

