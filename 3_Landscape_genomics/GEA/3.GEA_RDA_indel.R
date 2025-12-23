library(psych)    
library(vegan) 
library(adegenet)
library(LEA)
library(data.table)
library(parallel)
library(cols4all)

# ==============================================================================
# 1. 设置工作目录与数据加载
# ==============================================================================
setwd("/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/GEA_2025/RDA/indel")

# --- 加载基因型数据 ---
# 读取 .lfmm 文件
rawgeno <- fread("/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.missing0.05.plink.recodeA.lfmm")

# 读取 SNP ID
snp_ids <- fread("/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/indel/202_samples_indel_filtered.missing0.05.ID", 
                 header = FALSE)

# 检查维度
if (ncol(rawgeno) != nrow(snp_ids)) stop("SNP/INDEL 标识符的数量与基因型数据的列数不一致。")
colnames(rawgeno) <- snp_ids$V1

# --- 基因型插补 (Imputation) ---
# 使用众数(Mode)填补缺失值
gen.imp <- apply(rawgeno, 2, function(x) {
  if(any(is.na(x))) {
    modes <- as.numeric(names(which.max(table(x))))
    x[is.na(x)] <- modes
  }
  return(x)
})

# 保存插补后基因型数据
save(gen.imp, file = "gen.imp_indel.RData") 

# --- 加载环境数据 ---
env_data = read.csv("/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/GEA_2025/extracted_data/Climate_current_202samples.csv",
                    header = T, row.names = 1)
env_data = env_data[, 3:21] # 截取BIO1-19
str(env_data)

# 确保数据对齐 (假设顺序一致，若不一致请添加对齐代码)

# ==============================================================================
# 2. 构建 RDA 模型
# ==============================================================================

# 定义环境变量公式
env_formula <- gen.imp ~ BIO2 + BIO7 + BIO8 + BIO9 + BIO10 + BIO12 + BIO15 + BIO17 + BIO18

# 运行 RDA
of.rda <- rda(env_formula, data = env_data)
save(of.rda, file = "of.rda_indel.RData") 

# --- 关键步骤：共线性检查 (VIF) ---
## vif_res <- vif.cca(of.rda)
## print(vif_res)

# --- 关键步骤：R2 调整 ---
R2 <- RsquareAdj(of.rda)
print(R2)
## $r.squared
## [1] 0.1015557
## $adj.r.squared
## [1] 0.05944109

# 查看各轴的解释方差
eig_summary <- summary(eigenvals(of.rda, model = "constrained"))
eig_summary

# 保存特征值表
write.csv(eig_summary, file = "Importance_of_components_eigenvals_RDA_indel.csv")

# ==============================================================================
# 3. 显著性检验
# ==============================================================================

# 1. 检验全模型是否显著
## anova_global <- anova.cca(of.rda, parallel= 20) 
## print(anova_global)

# 2. 检验各个轴是否显著
## anova_axis <- anova.cca(of.rda, by="axis", parallel= 20)
## print(anova_axis)


# ==============================================================================
# 4. 离群位点筛选 (Outlier Detection)
# ==============================================================================

# 提取 SNP 在各轴上的载荷 (Loadings)
load.rda <- scores(of.rda, choices=c(1:3), display="species") 

# 定义筛选函数 (3倍标准差)
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# 对前三个轴分别筛选
cand1 <- outliers(load.rda[,1], 3)
cand2 <- outliers(load.rda[,2], 3)
cand3 <- outliers(load.rda[,3], 3)

# 合并所有候选位点名称 (去重)
names_list <- list(names(cand1), names(cand2), names(cand3))
combined_cand <- Reduce(union, names_list)

print(paste("共检测到候选 INDEL 数量:", length(combined_cand)))
## 共检测到候选 INDEL 数量: 28059

# --- 导出详细信息表格 ---
make_cand_df <- function(cand_vec, axis_num) {
  if(length(cand_vec) == 0) return(NULL)
  data.frame(axis = axis_num, snp = names(cand_vec), loading = unname(cand_vec))
}

df1 <- make_cand_df(cand1, 1)
df2 <- make_cand_df(cand2, 2)
df3 <- make_cand_df(cand3, 3)

cand_df_full <- rbind(df1, df2, df3)

# 保存 ID 列表
write.csv(data.frame(names = combined_cand), "RDA_cand_IDs_indel.csv", row.names = FALSE)

# 保存详细信息表
write.csv(cand_df_full, "RDA_cand_Details_indel.csv", row.names = FALSE)

# 检查是否有 INDEL 同时在多个轴上显著
duplicate_snps <- cand_df_full$snp[duplicated(cand_df_full$snp)]
if(length(duplicate_snps) > 0) {
  print(paste("有", length(unique(duplicate_snps)), "个 INDEL 在多个轴上均被识别为离群点"))
}
## 有 2054 个 INDEL 在多个轴上均被识别为离群点
