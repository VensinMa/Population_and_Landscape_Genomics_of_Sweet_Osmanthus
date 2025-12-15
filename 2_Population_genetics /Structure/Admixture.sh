# 1. 创建并进入结果目录
mkdir -p /home/vensin/workspace/snpcalling_wild/12.population_genetics/admixture/result
cd /home/vensin/workspace/snpcalling_wild/12.population_genetics/admixture/result

# 2.  并行运行 ADMIXTURE (从 K=20 到 K=2)
# -j 10: 同时运行 10 个 K 值任务
# --cv: 开启交叉验证 (用于通过 CV error 选择最佳 K)
# -s 123: 设置随机种子 (保证结果可复现)
# > log_K{}.out: 将包含 CV error 的日志重定向保存
seq 20 -1 2 | parallel -j 10 "admixture --cv -s 123 /home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.LD.pruned.bed {} > log_K{}.out" &

# 4. 确定最佳 K 值（CV error (Cross-Validation Error) 最小）
grep -h "CV" log_K*.out | sort -n -k 4

