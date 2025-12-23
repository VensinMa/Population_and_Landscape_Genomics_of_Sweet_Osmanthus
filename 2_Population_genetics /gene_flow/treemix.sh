## https://bitbucket.org/nygcresearch/treemix/downloads/

cd /home/vensin/software
wget https://bitbucket.org/nygcresearch/treemix/downloads/treemix-1.13.tar.gz
tar -zxvf treemix-1.13.tar.gz
cd treemix-1.13
./configure
make
make install

#########1.计算等位基因频率
# 1.1 vcf文件转换为tped格式，生成sample.tped和sample.tfam文件。
cd /home/vensin/workspace/snpcalling_wild/14.gene_flow/treemix
vcftools --gzvcf /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/208_samples_snp_filtered.LD.pruned.recode.vcf.gz --plink-tped --out 208samples

# 1.2 文件修改
# 208samples.tfam文件修改第一列数据为群体名称

# 1.3 提取文件sample.pop.cov，格式为：共三列，前两列与修改后的sample.tfam前两列一样，为群体ID和样本ID，第三列和第一列一致，tab分隔。
cat 208samples.tfam | awk '{print $1"\t"$2"\t"$1}' > 208samples.pop.cov

# 1.4 计算等位基因组的频率，生成plink.frq.strat和plink.nosex文件
plink --threads 4 --tfile 208samples --freq --allow-extra-chr --chr-set 23 --within 208samples.pop.cov

# 1.5 压缩等位基因频率文件
gzip plink.frq.strat

# 1.6 转换格式
#用treemix自带脚本进行格式转换，notes：输入输出都为压缩文件，将plink2treemix.py修改成与python3兼容
/home/vensin/anaconda3/envs/faststructure/bin/python /home/vensin/software/script/plink2treemix.py plink.frq.strat.gz 208samples.frq.gz

#########2.treemix推断基因流

#!/bin/bash

# ================= 配置区域 =================
INPUT="208samples.frq.gz"
OUTDIR="./root_MXL"
ROOT="O_MXL"
THREADS=16      # 设置同时运行的任务数 
# ===========================================

mkdir -p ${OUTDIR}

# 创建一个命令列表文件
CMD_FILE="treemix_commands.txt"
> ${CMD_FILE}

echo "正在生成任务列表..."

# 生成 0 到 15 个迁移事件，每个重复 100 次
for m in {0..15}; do
    for i in {1..100}; do
        # 生成随机种子，保证每次 bootstrap 都不一样
        SEED=$RANDOM
        
        # 构造输出前缀
        OUT_PREFIX="${OUTDIR}/mig${m}_iter${i}"
        
        # 构造完整命令
        # -k ${BLOCK_SIZE}: 设置 Block 大小
        # -s ${SEED}: 设置随机种子
        # > /dev/null 2>&1: 减少屏幕输出，防止刷屏
        echo "treemix -i ${INPUT} -m ${m} -bootstrap  -root ${ROOT} -s ${SEED} -o ${OUT_PREFIX} > /dev/null 2>&1" >> ${CMD_FILE}
    done
done

echo "任务列表生成完毕，共 $(wc -l < ${CMD_FILE}) 个任务。"
echo "开始并行执行，并发数: ${THREADS} ..."

# 使用 xargs 进行并行处理
cat ${CMD_FILE} | xargs -P ${THREADS} -I {} sh -c "{}"

echo "所有 TreeMix 任务执行完毕！"
# rm ${CMD_FILE} # 可选：删除临时任务文件
