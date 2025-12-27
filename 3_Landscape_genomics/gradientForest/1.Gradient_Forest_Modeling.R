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
# rm(list = ls()) # 清空环境变量
##############################  设置工作路径  #################################
workspace_dir <- "/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/gradientForest_2025/2.5m"
if (!dir.exists(workspace_dir)) dir.create(workspace_dir)
setwd(workspace_dir)
getwd() 
# 设置输入文件、输出结果和图片路径
result_dir <- "final_result"
picture_dir <- "final_picture"

if (!dir.exists(result_dir)) dir.create(result_dir)
if (!dir.exists(picture_dir)) dir.create(picture_dir)


#===========================    导入群体适应性位点等位基因频率数据
PopsMaf <- fread(file = "/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/gea_core_loci/GF_PopsMaf.csv")
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
AllEnvs = c("BIO2", "BIO7", "BIO8", "BIO9", "BIO10", 
            "BIO12", "BIO15", "BIO17", "BIO18")


# ====== 确保环境数据和基因数据的行名（群体）完全对齐 ======
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

# load(file = "gf.mod.10290.RData")
# 构建梯度森林模型，将环境数据和等位基因频率合并在一起
gf.mod <- gradientForest(Envs_Maf,
                         predictor.vars = colnames(CurrentEnvs),# 指定用作预测变量X的列名，即环境数据的列名
                         response.vars = colnames(PopsMaf),# 指定用作响应变量Y的列名，即等位基因频率数据的列名
                         ntree = 1000, maxLevel = maxLevel, trace = T,
                         corr.threshold = 0.50,  nbin = 1001, check.names = FALSE)

# 可以将保存 GF模型结果 gf.mod 保存到到一个文件 便于后续直接加载使用
save(gf.mod, file = "gf.mod.10290.RData")


########################################### 绘图 ###################################################
# 保存梯度森林模型的一些输出结果
# 包括响应变量 Y、预测变量 X、重要性值 imp.rsq、模型结果 result等
write.table(gf.mod$Y, file = paste0(result_dir, "/gf.mod.Y.txt"))
write.table(gf.mod$X, file = paste0(result_dir, "/gf.mod.X.txt"))
write.table(gf.mod$imp.rsq, file = paste0(result_dir, "/gf.mod.imp.rsq.txt"))
write.table(gf.mod$result, file = paste0(result_dir, "/gf.mod.result.txt"))
write.table(gf.mod$res.u, file = paste0(result_dir, "/gf.mod.res.u.txt"))
write.table(gf.mod$res, file = paste0(result_dir, "/gf.mod.res.txt"))

# 绘制预测变量重要性排序图
#生成空的PDF文件 #生成重要值排序 #保存生成的结果
pdf(file = paste0(picture_dir, "/gf.mod.CORE.Importance.pdf"), width = 8, height = 4) 
plot(gf.mod, plot.type = "Overall.Importance", 
     col = c(rep("grey",8), MaizePal::maize_pal("HighlandMAGIC", 3)),
     las = 2, cex.names = 0.8) 
dev.off() 


#======================= splits density plots 分割密度图plot.gradientForest ===
PredictEnvs = AllEnvs
PredictEnvs # [1] "BIO2"  "BIO7"  "BIO8"  "BIO9"  "BIO10" "BIO12" "BIO15" "BIO17" "BIO18"

## 可以先绘制单个环境因子的分割密度图，调整绘图参数
pdf(file = paste0(picture_dir, "/splits.density.plots_BIO2.pdf"), 
    width = 9, height = 6)
plot(gf.mod, plot.type= "S", imp.vars = "BIO2", leg.posn = "topleft", 
     cex.legend = 1, cex.axis = 1.2, cex.lab = 1.3, line.ylab = -1, 
     par.args = list(mgp = c(2, 0.5, 0), mar = c(3,3,1,1)))
dev.off()

#hist(gf.mod$X$BIO17, main="Histogram of BIO17", xlab="BIO17", breaks=20)

# 循环绘制 PredictEnvs 中所有环境变量的分割密度图
for (env_var in PredictEnvs) {
  pdf(file = paste0(picture_dir, "/GF_Splits_Density_", env_var, ".pdf"), 
      width = 9, height = 6)
  plot(gf.mod, plot.type = "S", imp.vars = env_var, leg.posn = "topright", 
       cex.legend = 1, cex.axis = 1.2, cex.lab = 1.3, line.ylab = -1, 
       par.args = list(mgp = c(2, 0.5, 0), mar = c(3, 3, 1, 1)))
  dev.off()
}

# 绘制 PredictEnvs 中所有环境变量的分割密度图到一个PDF
pdf(file = paste0(picture_dir, "/GF_Splits_Density_All_Variables.pdf"), width = 9, height = 6)
for (var in AllEnvsNames) {
  print(paste("正在绘制:", var))
  plot(gf.mod, plot.type = "S", imp.vars = var, 
       leg.posn = "topright", # 图例位置
       cex.legend = 1,      # 图例大小
       cex.axis = 1.2,          # 坐标轴字体大小
       cex.lab = 1.3,         # 标签字体大小
       line.ylab = -1,       # Y轴标签距离
       # 调整边距: 下、左、上、右
       par.args = list(mgp = c(2, 0.5, 0), mar = c(3, 2, 2, 1))) 
}
dev.off()

#======================= Cumulative Importance plots  ==========================

#显示整体（环境因子BIO1-19）组成的累积变化，其中变化发生在梯度上
pdf(file = paste0(picture_dir, "/GF_Cumulative_Importance_All_Variables_C.pdf"), 
    width = 6, height = 6)
plot(gf.mod, plot.type = "C", imp.vars = AllEnvsNames, show.species = F, 
     common.scale = F, cex.axis = 1, cex.lab = 1.2, line.ylab = 0.8, 
     par.args = list(mgp = c(1.7, 0.5, 0), mar = c(3, 2, 0.1, 0.5), 
                     omi = c(0, 0.3, 0.05, 0.1)))
dev.off()

# ## 单个环境因子
# # pdf(file = paste0(picture_dir, "/GF_Cumulative_Importance_BIO2.pdf"), width = 5, height = 7)
# plot(gf.mod, plot.type = "Cumulative.Importance", imp.vars = "BIO2", leg.nspecies = 5,
#      show.overall = T, legend = T, leg.posn = "topleft", common.scale = T, 
#      cex.lab = 1, cex.legend = 0.7, cex.axis = 1, line.ylab = 2,
#      par.args = list(mgp = c(2, 0.5, 0), mar = c(0.5, 0.2, 0.2, 1),
#                      omi = c(0.3, 0.7, 0.1, 0.1)))
# mgp 距轴的距离=c(轴标签，轴刻度标签，轴刻度)；内边距 mar=c(下，左，上，右)；外边距 omi
# dev.off()

# # 循环迭代每个环境因子
# for (env_factor in PredictEnvs) {
#   pdf_file <- paste0(picture_dir, "/GF_Cumulative_Importance_", env_factor, ".pdf")
#   pdf(file = pdf_file, width = 5, height = 7)
#   plot(gf.mod, plot.type = "Cumulative.Importance", imp.vars = env_factor, leg.nspecies = 5,
#        show.overall = T, legend = T, leg.posn = "topleft", common.scale = T, 
#        cex.lab = 1, cex.legend = 0.7, cex.axis = 1, line.ylab = 2,
#        par.args = list(mgp = c(2, 0.5, 0), mar = c(0.5, 0.2, 0.2, 1), 
#                        omi = c(0.3, 0.7, 0.1, 0.1)))
#   dev.off()
# }

# pdf(file = paste0(picture_dir, "/GF_Cumulative_Importance.pdf"), width = 8, height = 6)
# # plot.type = "C" 表示 Cumulative
# plot(gf.mod, plot.type = "C", imp.vars = AllEnvsNames, 
#      show.overall = TRUE, legend = TRUE, leg.posn = "topleft",
#      common.scale = TRUE, 
#      col = rainbow(length(AllEnvsNames)), # 或者使用你喜欢的配色
#      cex.lab = 1.2, cex.axis = 1.2, line.ylab = 0.9)
# dev.off()

#============================   R2 重要性 plots  ===============================

pdf(file = paste0(picture_dir, "/GF_R2_Performance_Histogram.pdf"), width = 8, height = 6)
# 提取 R2 数据,绘制直方图
r2_values <- gf.mod$result
hist(r2_values, breaks = 50, col = "steelblue", border = "white",
     main = "Distribution of SNP R-squared values",
     xlab = "R-squared (Predictability)", ylab = "Frequency (Number of SNPs)")
abline(v = mean(r2_values, na.rm = T), col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("Mean R2 =", round(mean(r2_values, na.rm = T), 4)), 
       col = "red", lty = 2, lwd = 2)
dev.off()

# 
pdf(file = paste0(picture_dir, "/GF_R2_Performance_rank.pdf"), width = 7, height = 7)
plot(gf.mod, plot.type = "Performance", show.names = T, horizontal = T, 
     cex.axis = 1.2, cex.labels = 0.8, line = 2, 
     par.args = list(mgp = c(0, 0.8, 0),
                     mar = c(4, 7, 2, 0.5), omi = c(0, 0, 0.1, 0.1)))
dev.off()

#============================   GF MAP—PCA  ===============================
# 从 CSV 文件中读取数据，该文件包含了研究区域的坐标点和current气候数据
All_current_xy_envs = fread("extracted_data/Climate_current_66239_grids.csv")
All_current_xy_envs = as.data.frame(All_current_xy_envs)

# 检查模型所需的变量是否都在All_current_xy_envs数据中
if(!all(PredictEnvs %in% colnames(All_current_xy_envs))) {
  stop("错误：网格数据缺少模型所需的某些环境变量！")
}

# 筛选列并移除缺失值
# 只保留 经纬度 + 预测变量
cols_needed <- c("lon", "lat", PredictEnvs)
All_current_xy_envs <- na.omit(All_current_xy_envs[, cols_needed])
cat(paste("最终用于预测的网格数量:", nrow(All_current_xy_envs), "\n"))

# predict() 会返回转换后的生物学重要性数值 (cumulative importance)
gf_trans <- predict(gf.mod, All_current_xy_envs[, PredictEnvs])

# 将坐标和转换后的数据合并
# 注意：这里 gf_trans 里的列名通常就是 PredictEnvs 的名字
Trns_grid <- cbind(All_current_xy_envs[, c("lon", "lat")], gf_trans)

# 再次检查是否有 NA (理论上上一步清理过，这里不应该有 NA，但为了保险)
if (any(is.na(Trns_grid))) {
  warning("检测到预测结果中存在 NA，正在移除...")
  Trns_grid <- na.omit(Trns_grid)
}

# 使用主成分分析（PCA）进行降维，选择基于gf重要值排序和排除自相关筛选出的PredictEnvs
# Trns_grid[, PredictEnvs]：只对“环境变量列”做 PCA，不能包含 lon 和 lat
# scale. = FALSE : 因为 GF 的输出已经是“重要性”量纲，如果 scale=TRUE 会导致不重要的变量噪声被放大
All_PCs <- prcomp(Trns_grid[, PredictEnvs], center = TRUE, scale. = FALSE)
# 输出 PCA 摘要
print(summary(All_PCs))

# 为地图绘制设置颜色，使用主成分的分数进行颜色映射
a1 <- All_PCs$x[, 1]
a2 <- All_PCs$x[, 2]
a3 <- All_PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
# 归一化颜色值到0-255范围
r <- (r - min(r)) / (max(r) - min(r)) * 255
g <- (g - min(g)) / (max(g) - min(g)) * 255
b <- (b - min(b)) / (max(b) - min(b)) * 255
# 创建包含坐标和颜色信息的数据框
grid <- All_current_xy_envs[, c("lon","lat")]
grid$R = r
grid$G = g
grid$B = b
# 获取主成分分析中环境因子的数量
nvs <- dim(All_PCs$rotation)[1]
nvs
vec <- PredictEnvs
lv <- length(vec)
lv
vind <- rownames(All_PCs$rotation) %in% vec
scal <- 60
xrng <- range(All_PCs$x[, 1], All_PCs$rotation[, 1]/scal) * 1.1
yrng <- range(All_PCs$x[, 2], All_PCs$rotation[, 2]/scal) * 1.1

############################  PC散点图  #################################
library(ggplot2)

# 准备主成分得分数据 (PC1 和 PC2)
pc_scores <- as.data.frame(All_PCs$x[, 1:2])
colnames(pc_scores) <- c("PC1", "PC2")

# 添加颜色数据 (RGB)
pc_scores$R <- r
pc_scores$G <- g
pc_scores$B <- b
pc_scores$Color <- rgb(pc_scores$R, pc_scores$G, pc_scores$B, maxColorValue = 255)

# 准备方向向量数据
arrows_data <- as.data.frame(All_PCs$rotation[, 1:2])
colnames(arrows_data) <- c("ArrowX", "ArrowY")
arrows_data$Var <- vec  # 添加环境因子名称
arrows_data[, c("ArrowX", "ArrowY")] <- arrows_data[, c("ArrowX", "ArrowY")] / scal # 对数据框的数值列进行缩放
jit = 0.002

# 提取 PC1 和 PC2 的解释度
summ <- summary(All_PCs)
# 提取第二行(Proportion of Variance)的第一列(PC1)和第二列(PC2)，并转换为百分比保留2位小数
pc1_var <- round(summ$importance[2, 1] * 100, 2)  # 对应你的 56.26
pc2_var <- round(summ$importance[2, 2] * 100, 2)  # 对应你的 30.59

pdf(file = paste0(picture_dir, "/GF_PCplot_PC1_PC2.pdf"), width = 7.5, height = 5.5)
ggplot() +
  # 绘制散点图
  geom_point(data = pc_scores, aes(x = PC1, y = PC2, color = Color), size = 1.2, shape = 16) +
  scale_color_identity() +  # 直接使用提供的颜色
  # 添加方向向量
  geom_segment(data = arrows_data, aes(x = 0, y = 0, xend = ArrowX, yend = ArrowY),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.4) +
  # 添加环境因子标签
  geom_text(data = arrows_data, aes(x = ArrowX + jit * sign(ArrowX), 
                                    y = ArrowY + jit * sign(ArrowY), 
                                    label = Var),
            size = 4, color = "black") +
  coord_fixed(ratio = 1) +  # 确保 x 和 y 轴比例一致
  theme_bw() + # 简洁的主题
  # 调整轴标签字体和颜色
  theme(
    axis.title.x = element_text(size = 16, color = "black"),  # x轴标签字体大小和颜色
    axis.title.y = element_text(size = 16, color = "black"),  # y轴标签字体大小和颜色
    axis.text = element_text(size = 14, color = "black"),     # x和y轴刻度文字大小
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次网格线
  ) + 
  labs(
    x = paste0("Principal Component 1 (", pc1_var, "%)"), 
    y = paste0("Principal Component 2 (", pc2_var, "%)"), 
    title = NULL
  )
dev.off()

############################  PC地图   #################################
library(ggplot2)

# 创建包含经纬度和颜色的绘图数据框
map_data <- as.data.frame(Trns_grid[, c("lon", "lat")])  # 提取经纬度数据
colnames(map_data) <- c("Longitude", "Latitude")
map_data$Color <- rgb(r, g, b, max = 255)  # 将颜色信息加入数据框

# 绘制主成分地图
pdf(file = paste0(picture_dir, "/GF_PCplot_MAP.pdf"), width = 7.5, height = 5.5)
ggplot(map_data, aes(x = Longitude, y = Latitude)) +
  # 绘制点图
  geom_point(aes(color = Color), shape = 17,size = 0.2) +  # size 控制点的大小
  scale_color_identity() +  # 直接使用 RGB 颜色
  coord_fixed(ratio = 1) +  # 确保经纬度比例一致
  scale_y_continuous(limits = c(19, 35)) +  # 设置 Y 轴范围，
  theme_bw() +  # 使用简洁的主题
  # 设置标题和轴标签
  labs(
    title = NULL,  # 添加主标题
    x = "Longitude", 
    y = "Latitude") +
  # 调整主题样式
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # 标题居中加粗
    axis.title = element_text(size = 16, color = "black"),  # 坐标轴 标签字体大小颜色
    axis.text = element_text(size = 14, color = "black"),  # 坐标轴 刻度字体大小颜色
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    panel.background = element_blank())  # 移除背景填充
dev.off()

##########################  保存PC颜色信息 ####################################

# 将RGB颜色转换为ArcGIS使用的格式
greencols=rgb(r,g,b, max=255)
greencols2=col2rgb(greencols)
greencols3=t(greencols2)

# 创建包含颜色信息的数据集
gradients=cbind(Trns_grid[, 1:2],greencols3)
gradients$color=greencols

# 将颜色数据保存为CSV文件
write.csv(gradients,file = paste0(result_dir, "/all_gradients4arcgis.csv"),
          row.names=F,quote=F)

#################################  THE END  #####################################



