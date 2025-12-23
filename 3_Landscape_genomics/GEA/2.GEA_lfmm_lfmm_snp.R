# ==============================================================================
# 1. 环境设置与包安装
# ==============================================================================

if(!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if(!requireNamespace("lfmm", quietly = TRUE))    remotes::install_github("bcm-uga/lfmm")
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("qvalue", quietly = TRUE)) BiocManager::install("qvalue")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")

library(lfmm)
library(qvalue)  
library(data.table)

# 设置工作目录 (请确保路径正确)
setwd("/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/GEA_2025/lfmm/snp")
cat("当前工作目录:", getwd(), "\n")

# ==============================================================================
# 2. 数据读取与预处理
# ==============================================================================

# --- 读取环境数据 (解释变量 X) ---
X_df = read.csv("../../extracted_data/Climate_current_202samples.csv", 
                header = T, row.names = 1)
X_df = X_df[, 3:21] # 截取 BIO1 - BIO19

# 获取环境因子名称 (用于后续给结果列命名)
env_names <- colnames(X_df)

# 关键：转换为矩阵 (Matrix)，lfmm 计算需要矩阵格式
X <- as.matrix(X_df)
cat("环境数据维度:", dim(X), "\n")
## 环境数据维度: 202 19 

# --- 读取基因型数据 (响应变量 Y) ---
# 读取 .lfmm 文件 (没有表头)
Y_path <- "/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.missing0.05.plink.recodeA.lfmm"
Y <- fread(Y_path, header = FALSE)

# 关键：转换为矩阵
Y <- as.matrix(Y)
cat("基因型数据维度:", dim(Y), "\n")
## 基因型数据维度: 202 10590591 

# 安全检查：确保 X 和 Y 的样本量一致
if(nrow(X) != nrow(Y)) stop("错误：环境数据 (X) 与 基因型数据 (Y) 的样本数(行数)不匹配！")

# --- 读取 snp ID ---
snp_ids <- fread("/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202_samples_snp_filtered.missing0.05.ID", 
                 header = FALSE)

# 安全检查：确保 snp 数量 与 ID 数量一致
if (ncol(Y) != nrow(snp_ids)) stop("错误：Y 矩阵的列数与 snp ID 文件的行数不匹配！")

# ==============================================================================
# 3. 运行 LFMM 模型
# ==============================================================================

K = 4  # 最佳分群数 

cat("正在运行 LFMM Ridge (K =", K, ")...\n")
mod.lfmm <- lfmm::lfmm_ridge(Y = Y, 
                             X = X, 
                             K = K)
saveRDS(mod.lfmm, file = "mod.lfmm.rds")
# mod.lfmm <- readRDS("mod.lfmm.rds") 

cat("正在进行统计检验 (GIF 校正)...\n")
pv <- lfmm::lfmm_test(Y = Y, X = X, 
                      lfmm = mod.lfmm, 
                      calibrate = "gif")
saveRDS(pv, file = "pv.rds")
# pv <- readRDS("pv.rds") 

# ==============================================================================
############# 如内存不足，可删除部分环境变量，释放内存 #########################\
# 128G内存在这里占用70%左右
rm(mod.lfmm)  # 删除模型对象
rm(Y)         # 删除基因型大矩阵
gc()          # 强制进行垃圾回收，释放内存回系统
print(ls())   # 释放后内存占用40% 左右
# ==============================================================================
# 4. 提取与保存 P 值 (P-values)
# ==============================================================================

# 提取 P 值
raw_pvalues <- pv$pvalue
calibrated_pvalues <- pv$calibrated.pvalue

# 查看基因组膨胀因子 (GIF)
cat("基因组膨胀因子 (GIF):\n")
print(pv$gif)

# 关键优化：给矩阵赋予列名 (环境因子名称)，否则结果只有数字，分不清是哪个环境因子
colnames(raw_pvalues) <- env_names
colnames(calibrated_pvalues) <- env_names

# 准备数据框
# 注意：data.table::fwrite 对 row.names=TRUE 支持一般，建议显式创建 SNP_ID 列
df_raw <- data.frame(SNP_ID = snp_ids$V1, raw_pvalues)
df_cal <- data.frame(SNP_ID = snp_ids$V1, calibrated_pvalues)

# 保存合并后的数据
fwrite(df_raw, file = "snp_raw_pvalues_merged.csv", sep = ",")
fwrite(df_cal, file = "snp_calibrated_pvalues_merged.csv", sep = ",")

# ==============================================================================
# 5. Q 值校正 (FDR) 与最终保存
# ==============================================================================
cat("正在计算 Q 值 (FDR)...\n")

# 初始化矩阵存储所有 Q 值
all_qvalues <- matrix(NA, nrow = nrow(calibrated_pvalues), ncol = ncol(calibrated_pvalues))
colnames(all_qvalues) <- env_names

# 循环处理每个环境因子
for (i in 1:ncol(calibrated_pvalues)) {
  
  # 提取当前环境因子的 P 值
  pvals <- calibrated_pvalues[, i]
  
  # 使用 try() 防止某个因子计算出错导致循环中断
  try({
    qobj <- qvalue(p = pvals)
    all_qvalues[, i] <- qobj$qvalues
    
    # 如需要每个因子单独存q值文件，可取消下面注释
    tmp_df <- data.frame(SNP_ID = snp_ids$V1, qvalues = qobj$qvalues)
    fwrite(tmp_df, file = paste0("qvalues_env_", env_names[i], ".csv"))
  })
}

# 合并 SNP ID 和 Q 值
df_final_q <- data.frame(SNP_ID = snp_ids$V1, all_qvalues)

# 保存最终汇总文件
fwrite(df_final_q, file = "all_EVs_qvalues_combined.csv", sep = ",")

cat("分析完成，所有文件已保存。\n")

# ==============================================================================
# 6. 显著性位点筛选 (FDR < 0.05 和 FDR < 0.01)
# ==============================================================================

# 确保列名存在 (承接上一步代码)
if(is.null(colnames(calibrated_pvalues))) colnames(calibrated_pvalues) <- env_names
env_factor_names <- colnames(calibrated_pvalues)
num_env_factors <- ncol(calibrated_pvalues)

# 初始化列表用于存储结果
threshold_list <- list()
significant_counts_list <- list() # 使用 list 动态存储，最后合并
all_sig_snps_vector <- c() # 用于存储所有显著位点的向量 (用于求并集)

cat("正在筛选显著位点并统计...\n")

for (i in 1:num_env_factors) {
  
  # 获取当前环境因子的名称和 P 值
  curr_env <- env_factor_names[i]
  pvals <- calibrated_pvalues[, i]
  
  # --- 1. 计算 Q 值 ---
  # 注意：如果在上一步代码中已经计算过 all_qvalues，直接调用即可，无需重复计算
  # 这里为了代码独立性保留计算过程，但加入了 tryCatch
  qvals <- tryCatch({
    qvalue(p = pvals)$qvalues
  }, error = function(e) {
    cat("Warning: qvalue calculation failed for", curr_env, "\n")
    return(rep(NA, length(pvals)))
  })
  
  if(all(is.na(qvals))) next # 如果计算失败则跳过该循环
  
  # --- 2. 寻找 P 值阈值 (更严谨的方法) ---
  # 方法：找到所有 q < 0.05 的位点中，P 值最大的那个，作为该环境因子的阈值
  
  # 针对 FDR < 0.05
  sig_idx_05 <- which(qvals < 0.05)
  if(length(sig_idx_05) > 0) {
    threshold_p_05 <- max(pvals[sig_idx_05])
    actual_q_05 <- max(qvals[sig_idx_05])
  } else {
    threshold_p_05 <- NA
    actual_q_05 <- NA
  }
  
  # 针对 FDR < 0.01
  sig_idx_01 <- which(qvals < 0.01)
  if(length(sig_idx_01) > 0) {
    threshold_p_01 <- max(pvals[sig_idx_01])
    actual_q_01 <- max(qvals[sig_idx_01])
  } else {
    threshold_p_01 <- NA
    actual_q_01 <- NA
  }
  
  # 保存阈值信息
  threshold_list[[i]] <- data.frame(
    Environment_Factor = curr_env,
    Threshold_P_05 = threshold_p_05,
    Max_Q_below_05 = actual_q_05,
    Threshold_P_01 = threshold_p_01,
    Max_Q_below_01 = actual_q_01
  )
  
  # --- 3. 保存该环境因子的所有结果 ---
  # 包含 SNP ID, P值, Q值
  res_df <- data.frame(SNPid = snp_ids$V1, pvalues = pvals, qvalues = qvals)
  # 建议：仅保存显著的或者保存一个汇总大文件，避免生成几十个 csv 小文件。
  # 如果确实需要每个因子一个文件：
  fwrite(res_df, file = paste0("all_results_", curr_env, ".csv"), row.names = FALSE)
  
  # --- 4. 提取显著性位点 (FDR < 0.05) ---
  if (length(sig_idx_05) > 0) {
    
    # 提取信息
    sig_snps <- snp_ids$V1[sig_idx_05]
    sig_p <- pvals[sig_idx_05]
    sig_q <- qvals[sig_idx_05]
    
    # 保存显著结果详情
    sig_df <- data.frame(SNPid = sig_snps, pvalues = sig_p, qvalues = sig_q)
    fwrite(sig_df, file = paste0("significant_results_", curr_env, ".csv"), row.names = FALSE)
    
    # 仅保存 SNP ID 列表 (用于后续提取基因等)
    fwrite(data.frame(SNPid = sig_snps), file = paste0("significant_snp_ids_", curr_env, ".csv"), row.names = FALSE, col.names = FALSE)
    
    # 添加到总向量中 (用于后续求并集)
    all_sig_snps_vector <- c(all_sig_snps_vector, sig_snps)
    
    # 记录数量
    count_val <- length(sig_snps)
    
  } else {
    count_val <- 0
  }
  
  significant_counts_list[[i]] <- data.frame(Environment_Factor = curr_env, Significant_SNP_Count = count_val)
}

# --- 5. 汇总与保存统计信息 ---

# 合并并保存阈值表
threshold_df <- do.call(rbind, threshold_list)
fwrite(threshold_df, file = "environment_factors_thresholds_summary.csv")

# 合并并保存显著位点计数表
counts_df <- do.call(rbind, significant_counts_list)
fwrite(counts_df, file = "significant_snp_counts.csv")

# ==============================================================================
# 7. 统计所有环境因子关联位点的并集 (Union)
# ==============================================================================

# 去除重复，得到唯一的显著位点列表
union_sig_snps <- unique(all_sig_snps_vector)
union_count <- length(union_sig_snps)

# 保存并集结果
if (union_count > 0) {
  union_df <- data.frame(SNPid = union_sig_snps)
  fwrite(union_df, file = "LFMM_significant_snp_ids_union.csv", row.names = FALSE)
  cat("\n所有环境因子下显著位点(FDR < 0.05)的并集数量为:", union_count, "\n")
} else {
  cat("\n没有检测到任何显著关联位点 (FDR < 0.05)。\n")
}

## 所有环境因子下显著位点(FDR < 0.05)的并集数量为: 346683 
cat("筛选完成。\n")
