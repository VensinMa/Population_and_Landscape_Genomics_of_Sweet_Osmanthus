# 1. 创建并进入结果目录
mkdir -p /home/vensin/workspace/snpcalling_wild/12.population_genetics/admixture/result
cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/admixture/result

# 2.  并行运行 ADMIXTURE (从 K=20 到 K=2)
# -j 10: 同时运行 10 个 K 值任务
# --cv: 开启交叉验证 (用于通过 CV error 选择最佳 K)
# -s 123: 设置随机种子 (保证结果可复现)
# > log_K{}.out: 将包含 CV error 的日志重定向保存
seq 2 20 | parallel -j 10 "admixture --cv -s 123 /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.bed {} > log_K{}.out" &

# 4. 确定最佳 K 值（CV error (Cross-Validation Error) 最小）
grep -h "CV" log_K*.out | sort -n -k 4
## (base) vensin@Ubuntu24:~/workspace/snpcalling_wild/12.population_genetics/admixture/result$ grep -h "CV" log_K*.out | sort -n -k 4
## CV error (K=6): 0.43273
## CV error (K=5): 0.43483
## CV error (K=4): 0.43840
## CV error (K=3): 0.44112
## CV error (K=7): 0.44222
## CV error (K=2): 0.44499
## CV error (K=8): 0.44512
## CV error (K=9): 0.44775
## CV error (K=10): 0.45419
## CV error (K=11): 0.45523
## CV error (K=12): 0.45850
## CV error (K=13): 0.46217
## CV error (K=14): 0.48551
## CV error (K=15): 0.48637
## CV error (K=16): 0.48728
## CV error (K=17): 0.50524
## CV error (K=19): 0.51016
## CV error (K=20): 0.52380
## CV error (K=18): 0.52457

