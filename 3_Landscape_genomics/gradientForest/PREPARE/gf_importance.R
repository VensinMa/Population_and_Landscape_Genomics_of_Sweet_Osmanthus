# 加载所需的R包
# install.packages("conformal", repos="http://R-Forge.R-project.org")
# install.packages("extendedForest", repos="http://R-Forge.R-project.org")
# install.packages("gradientForest", repos="http://R-Forge.R-project.org")
# devtools::install_github("AndiKur4/MaizePal")
### sudo apt-get install -y libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
### sudo apt-get install -y cmake libudunits2-dev libgdal-dev libgeos-dev libproj-dev libabsl-dev


sink("my_log_20251220.txt")
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

#############################  GradientForest 模型的构建 #########################################

# 导入群体适应性位点等位基因频率数据
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

# 构建梯度森林模型
# 计算maxLevel，用于树的最大深度
maxLevel <- log2(0.368 * nrow(PopsMaf) / 2)
Envs_Maf = cbind(CurrentEnvs[, AllEnvsNames], PopsMaf)

# 构建梯度森林模型，将环境数据和等位基因频率合并在一起
gf.mod <- gradientForest(Envs_Maf,
                         predictor.vars = colnames(CurrentEnvs),# 指定用作预测变量X的列名，即环境数据的列名
                         response.vars = colnames(PopsMaf),# 指定用作响应变量Y的列名，即等位基因频率数据的列名
                         ntree = 1000, maxLevel = maxLevel, trace = T,
                         corr.threshold = 0.50,  nbin = 1001, check.names = FALSE)

# 可以将保存 GF模型结果 gf.mod 保存到到一个文件 便于后续直接加载使用
save(gf.mod, file = "gf.mod.165692.RData")
# load(file = "gf.mod.165692.RData")

########################################### 绘图 ###################################################

# 绘制预测变量重要性排序图
#生成空的PDF文件 #生成重要值排序 #保存生成的结果
pdf(file = paste0(picture_dir, "/gf.mod.165692.Importance.pdf"), width = 8, height = 4) 
plot(gf.mod, plot.type = "Overall.165692.Importance", 
     col = c(rep("grey",8), MaizePal::maize_pal("HighlandMAGIC", 3)),
     las = 2, cex.names = 0.8) 
dev.off() 

# 提取重要性
importance_table <- data.frame(Variable = names(gf.mod$overall.imp), 
                               Importance = gf.mod$overall.imp)
importance_table <- importance_table[order(importance_table$Importance, decreasing = F), ]

# 设置因子水平以便画图排序
importance_table$Variable <- factor(importance_table$Variable, levels = importance_table$Variable)

# 绘制重要性条形图
p_imp <- ggplot(importance_table, aes(x = Importance, y = Variable)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Environmental Variable Importance", x = "Weighted Importance (R2)", y = "Variable")

ggsave(filename = file.path(picture_dir, "GF_Variable_Importance.pdf"), plot = p_imp, width = 6, height = 8)
