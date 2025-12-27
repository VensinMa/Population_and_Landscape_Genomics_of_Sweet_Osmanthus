# ==============================================================================
# 安装加载所需的R包，设置工作路径
# ==============================================================================

# install.packages("conformal", repos="http://R-Forge.R-project.org")
# install.packages("extendedForest", repos="http://R-Forge.R-project.org")
# install.packages("gradientForest", repos="http://R-Forge.R-project.org")
# devtools::install_github("AndiKur4/MaizePal")
### sudo apt-get install -y libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
### sudo apt-get install -y cmake libudunits2-dev libgdal-dev libgeos-dev libproj-dev libabsl-dev


library(gradientForest)
library(MaizePal)
library(data.table)
library(gdm)
library(dplyr)
library(tidyverse)
library(raster)
library(fields)
library(geosphere)
library(gdm)
library(foreach)
library(parallel)
library(doParallel)
library(gradientForest)
library(fields)
library(ggplot2)
library(sf)

getwd()
setwd("/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/gradientForest_2025")
getwd() 

# 设置输入文件、输出结果和图片路径
result_dir <- "result"
picture_dir <- "picture"

if (!dir.exists(result_dir)) dir.create(result_dir)
if (!dir.exists(picture_dir)) dir.create(picture_dir)

# ==============================================================================
# 导入等位基因频率数据与环境数据
# ==============================================================================

PopsMaf <- fread(file = "input/GF_PopsMaf.csv")
setnames(PopsMaf, old = "V1", new = "pop") # 第一列是群体名称但没有列名 将默认列名V1更改为pop
setDF(PopsMaf) # 将data.table转换为data.frame，因为行名是data.frame的特性
rownames(PopsMaf) <- PopsMaf[[1]] # 将第一列作为行名
PopsMaf <- PopsMaf[,-1] # 移除已经设置为行名的列
PopsMaf <- PopsMaf[order(rownames(PopsMaf)), ]
PopsMaf <- PopsMaf %>%
  dplyr::select(where(~ !any(is.na(.))))  # 去除缺失位点
dim(PopsMaf)

# df <- read.csv("/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/GEA_2025/extracted_data/Climate_current_202samples.csv")
# df_pop <- df %>%
#   mutate(Population = sub("_[^_]+$", "", ind)) %>%
#   group_by(Population) %>%
#   summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
# print(head(df_pop))
# write.csv(df_pop, "input/32pop_means_env_vars.csv", row.names = FALSE)
# ==============================================================================

# 读取群体当前环境数据
CurrentEnvs <- read.csv('input/32pop_means_env_vars.csv', 
                        header = T, row.names = 1)

# 选择最终保留的用于构建GF模型的预测环境因子
AllEnvs = c("BIO1", "BIO2", "BIO3", "BIO4", "BIO5", "BIO6", "BIO7", "BIO8", 
            "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", 
            "BIO16", "BIO17", "BIO18", "BIO19")


# === 【关键步骤】确保环境数据和基因数据的行名（群体）完全对齐 ===
common_pops <- intersect(rownames(PopsMaf), rownames(CurrentEnvs))
PopsMaf <- PopsMaf[common_pops, ]
CurrentEnvs <- CurrentEnvs[common_pops, AllEnvs] # 同时只保留需要的环境因子
rownames(PopsMaf)

if (!all(rownames(PopsMaf) == rownames(CurrentEnvs))) {
  stop("错误：环境数据和基因数据的群体顺序不一致！请检查行名。")
} else {
  print("检测通过：基因数据与环境数据群体顺序一致。")
}

print(paste("最终用于建模的群体数量:", length(common_pops)))
print(paste("最终用于建模的 SNP 数量:", ncol(PopsMaf)))
CurrentEnvs <- CurrentEnvs[, AllEnvs]
AllEnvsNames = colnames(CurrentEnvs)
print(AllEnvsNames)

# ==============================================================================
# GF模型构建：批量随机抽样构建 GF 模型 (10次循环)
# ==============================================================================

# 设置参数
n_repeats <- 10       # 循环次数
n_snps_sample <- 10000 # 每次抽取的 SNP 数量
out_path <- paste0(picture_dir, "/Bootstrap_10k_SNPs") # 结果保存子目录

if(!dir.exists(out_path)) dir.create(out_path)

# 创建一个列表来存储每次循环的变量重要性数值（为了后续画汇总图）
importance_list <- list()

# 计算树的最大深度 (这个只与群体数量有关，循环中是不变的)
maxLevel <- log2(0.368 * nrow(PopsMaf) / 2)

print(paste("开始运行随机抽样验证，共", n_repeats, "次循环..."))

for (i in 1:n_repeats) {
  
  # 1. 设置随机种子
  set.seed(2025 + i) 
  
  # 2. 随机抽取 SNP
  if(ncol(PopsMaf) < n_snps_sample) {
    warning("总 SNP 数量少于 10000，将使用全部 SNP。")
    sample_idx <- 1:ncol(PopsMaf)
  } else {
    sample_idx <- sample(ncol(PopsMaf), n_snps_sample)
  }
  
  # 3. 构建本次循环的输入数据
  Current_Maf_Sub <- PopsMaf[, sample_idx]
  Envs_Maf_Sub <- cbind(CurrentEnvs[, AllEnvsNames], Current_Maf_Sub)
  
  print(paste0(">>> 正在运行第 ", i, "/", n_repeats, " 次模型 (SNPs: ", length(sample_idx), ")..."))
  
  # 4. 构建 Gradient Forest 模型
  gf.mod.sub <- gradientForest(Envs_Maf_Sub,
                               predictor.vars = AllEnvsNames,
                               response.vars = colnames(Current_Maf_Sub),
                               ntree = 500, 
                               maxLevel = maxLevel, 
                               trace = F, 
                               corr.threshold = 0.50, 
                               nbin = 1001, 
                               check.names = FALSE)
  
  # 5. 提取并保存重要性数据
  imp_df <- data.frame(Variable = names(gf.mod.sub$overall.imp),
                       Importance = gf.mod.sub$overall.imp,
                       Run = paste0("Run_", i))
  importance_list[[i]] <- imp_df
  
  # 6. 绘制并保存单次循环的重要性图
  pdf(file = paste0(out_path, "/Importance_Run_", i, ".pdf"), width = 8, height = 6)
  
  # ---【修改点】删除 main 参数，避免冲突 ---
  plot(gf.mod.sub, 
       plot.type = "Overall.Importance",
       # main = ...,  <-- 这一行删掉了
       col = "steelblue", 
       las = 2, 
       cex.names = 0.8)
  
  # 如果你非常想要标题，可以在 plot 之后用 title() 函数追加（但可能会重叠，不建议）
  # title(main = paste0("Run ", i)) 
  
  dev.off()
  
  # 显式清理内存
  gc()
}

print("10 次循环运行完毕！")

save(importance_list, file = "importance_list.RData")

# ==============================================================================
# 汇总可视化：绘制变量重要性的箱线图
# ==============================================================================

# load(file = "importance_list.RData")
# 1. 将列表合并为一个大数据框
library(ggplot2)
all_importance <- do.call(rbind, importance_list)

# 2. 计算每个变量的平均重要性，用于排序
var_order <- all_importance %>%
  group_by(Variable) %>%
  summarise(Mean_Imp = mean(Importance)) %>%
  arrange(Mean_Imp) %>% # 升序排列，ggplot flip后变降序
  pull(Variable)

# 3. 将 Variable 转换为因子，指定排序
all_importance$Variable <- factor(all_importance$Variable, levels = var_order)

# 4. 使用 ggplot2 绘图
p_box <- ggplot(all_importance, aes(x = Variable, y = Importance)) +
  # --- 【新增】添加 T 型误差棒 ---
  # width 控制横杠的宽度，放在 boxplot 之前以保证图层层级正确
  stat_boxplot(geom = "errorbar", width = 0.3) + 
  
  # --- 原有的箱线图 ---
  geom_boxplot(fill = "lightblue", outlier.shape = NA, alpha = 0.7) + 
  
  # --- 原有的散点 ---
  geom_jitter(width = 0.2, size = 1, color = "darkblue", alpha = 0.6) + 
  
  # --- 坐标轴翻转与美化 ---
  coord_flip() + 
  theme_bw() +
  labs(title = "Stability of Predictor Importance (10 Runs x 10k SNPs)",
       x = "Environmental Variables",
       y = "Weighted Importance (R2)") +
  theme(axis.text.y = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5)) # 标题居中

# 显示图片
print(p_box)

# 保存图片
ggsave(filename = paste0(out_path, "/Summary_Importance_Boxplot_with_Errorbars.pdf"), 
       plot = p_box, width = 8, height = 6)




# ==============================================================================
# 0. 准备工作
# ==============================================================================
if(!require(corrplot)) install.packages("corrplot")
library(corrplot)
library(psych) 

# 读取数据
env_data <- read.csv("/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/GEA_2025/extracted_data/Climate_current_202samples.csv", header = TRUE, row.names = 1)

# 提取环境数据列 (假设在第3到21列)
raw_env_vars <- env_data[, 3:21]

# ---【关键步骤】强制按 BIO1 到 BIO19 的数值顺序排列列名 ---
desired_order <- paste0("BIO", 1:19) 
valid_cols <- intersect(desired_order, colnames(raw_env_vars))
raw_env_vars <- raw_env_vars[, valid_cols]

print("已将变量按 BIO1 - BIO19 顺序排列。")

# ==============================================================================
# 2. 计算皮尔逊相关系数矩阵
# ==============================================================================

# 直接使用 raw_env_vars 计算相关性
cor_matrix <- cor(raw_env_vars, method = "pearson", use = "pairwise.complete.obs")

# ==============================================================================
# 3. 绘图 (Heatmap)
# ==============================================================================

# 设置 PDF 输出
pdf(file = "picture/Environment_Correlation_Plot_Ordered.pdf", width = 10, height = 10)

# 绘制热图
corrplot(cor_matrix, 
         method = "color",       # 使用颜色填充
         type = "upper",         # 只显示上半部分
         
         order = "original",     # 使用原始顺序
         
         addCoef.col = "black",  # 显示相关系数数值
         
         # ---【修改点 1：设置显示3位小数】---
         number.digits = 3,      
         # --------------------------------
         
         tl.col = "black",       # 标签颜色
         tl.srt = 45,            # 标签旋转45度
         
         # 如果显示3位小数后数字太挤，可以适当调小字体，比如改成 0.5
         number.cex = 0.4,       
         
         diag = FALSE)           # 不显示对角线

dev.off()
print("相关性热图已保存至 picture/Environment_Correlation_Plot_Ordered.pdf")


# ==============================================================================
# 4. 筛选高相关性变量 (打印检查)
# ==============================================================================

# 获取相关性绝对值大于 0.8 的变量对
high_cor <- which(abs(cor_matrix) > 0.8, arr.ind = TRUE)

# 过滤掉对角线元素以及重复的组合
high_cor <- high_cor[high_cor[,1] < high_cor[,2], ]

# 打印结果
cat("\n=== 相关性绝对值 > 0.8 的变量对 ===\n")
if (nrow(high_cor) > 0) {
  for (i in 1:nrow(high_cor)) {
    row_idx <- high_cor[i, 1]
    col_idx <- high_cor[i, 2]
    
    var1 <- rownames(cor_matrix)[row_idx]
    var2 <- colnames(cor_matrix)[col_idx]
    cor_value <- cor_matrix[row_idx, col_idx]
    
    # ---【修改点 2：控制台输出也保留3位小数 (%.3f)】---
    cat(sprintf("Variables: %-6s and %-6s | Correlation: %.3f\n", var1, var2, cor_value))
  }
} else {
  cat("No variable pairs have a correlation above 0.8.\n")
}
