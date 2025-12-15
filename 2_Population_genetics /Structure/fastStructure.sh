# 1. 进入目录
mkdir -p /home/vensin/workspace/snpcalling_wild/12.population_genetics/faststructure/result
cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/faststructure

# 2. 运行 fastStructure (从 K=20 到 K=2)
# 使用 seq 20 -1 2
seq 2 20 | parallel -j 10 "structure.py -K {} --input=/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned --output=/home/vensin/workspace/snpcalling_wild/12.population_genetics/faststructure/result/LD_faststructure_K_{} --cv=5 --prior=simple --seed=123  > /home/vensin/workspace/snpcalling_wild/12.population_genetics/faststructure/result/log_K{}.log 2>&1" 


# 等上面所有任务跑完后运行    
cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/faststructure/result
chooseK.py --input=/home/vensin/workspace/snpcalling_wild/12.population_genetics/faststructure/result/LD_faststructure_K
## Model complexity that maximizes marginal likelihood = 4
## Model components used to explain structure in data = 7

