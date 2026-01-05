## cd /home/vensin/software
## conda create -n easySFS && conda activate easySFS
## git clone https://github.com/isaacovercast/easySFS.git

######################################################################## 
mkdir -p  /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/easySFS/40samples_folded
cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/easySFS/40samples_folded
easySFS.py -i  /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.recode.vcf -p ../40samples.pop -a  --preview

# --proj 顺序必须严格对应: Processing 4 populations - ['Southwest-Yunnan', 'Southwest-Guizhou', 'Central', 'East']
easySFS.py -i /home/vensin/workspace/snpcalling_wild/13.genetic_load/ancestral_allele/205_samples_snp_filtered.nomissing.recode.vcf \
  -p ../40samples.pop \
  -a \
  --proj 20,20,20,20 \
  -o output_folded_40samples

######################################################################## 
mkdir -p  /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/Stairwayplot2/40samples_folded_Intergenic_5kb
cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/Stairwayplot2/40samples_folded_Intergenic_5kb
easySFS.py -i  /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/Stairwayplot2/40samples_folded_Intergenic_5kb/205_samples_neutral_Intergenic_5kb.recode.vcf \
      -p /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/easySFS/40samples.pop -a  --preview

# --proj 顺序必须严格对应: Processing 4 populations - ['Southwest-Yunnan', 'Southwest-Guizhou', 'Central', 'East']
easySFS.py -i /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/Stairwayplot2/40samples_folded_Intergenic_5kb/205_samples_neutral_Intergenic_5kb.recode.vcf \
  -p /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/easySFS/40samples.pop \
  -a \
  --proj 20,20,20,20 \
  -o output_folded_40samples_Intergenic_5kb

######################################################################## 
mkdir -p  /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/Stairwayplot2/40samples_folded_Intergenic_5kb_re_sample

cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/Stairwayplot2/40samples_folded_Intergenic_5kb_re_sample
easySFS.py -i  /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/Stairwayplot2/40samples_folded_Intergenic_5kb/205_samples_neutral_Intergenic_5kb.recode.vcf \
      -p /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/Stairwayplot2/40samples_folded_Intergenic_5kb_re_sample/40samples.re.pop -a  --preview

# --proj 顺序必须严格对应: Processing 4 populations - ['Southwest-Yunnan', 'Southwest-Guizhou', 'Central', 'East']
easySFS.py -i /home/vensin/workspace/snpcalling_wild/12.population_genetics/Demographic_History/Stairwayplot2/40samples_folded_Intergenic_5kb/205_samples_neutral_Intergenic_5kb.recode.vcf \
  -p ./40samples.re.pop \
  -a \
  --proj 20,20,20,20 \
  -o output_folded_40samples.re_Intergenic_5kb



# 1. 设定 classpath (指向 stairway_plot_es 文件夹)
CP="/home/vensin/software/stairway_plot_v2.2/stairway_plot_es"

# 2. 批量生成 .sh 运行脚本
for file in *.blueprint; do
    echo "Processing $file ..."
    java -cp "$CP" Stairbuilder "$file"
done

# 3. 批量执行计算脚本
for script in *.blueprint.sh; do
    echo "Starting analysis for $script ..."
    # bash $script
    # 使用 nohup 后台运行，防止断线中断，同时记录日志
    nohup bash "$script" > "${script}.log" 2>&1 &
done
