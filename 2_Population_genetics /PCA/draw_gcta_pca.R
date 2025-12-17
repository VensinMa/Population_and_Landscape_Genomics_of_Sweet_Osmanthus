# 1. 安装及加载所需R包
packages <- c("ggplot2", "ggrepel", "ggsci", "data.table", "dplyr", "tidyr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
library(ggplot2)
library(ggrepel)
library(ggsci)
library(data.table)
library(dplyr)
library(tidyr)

# 2. 输入文件路径设置
eigen_vec_file <- "202_samples_gcta_pca_result.eigenvec" # PCA 坐标文件
eigen_val_file <- "202_samples_gcta_pca_result.eigenval" # PCA 特征值文件 (用于计算解释率)
pop_file       <- "202samples.pop"                       # 分组文件

# 输出设置
outpre <- "PCA_out_figure"
x_pc <- 1 # X轴使用第几主成分 (PC1)
y_pc <- 2 # Y轴使用第几主成分 (PC2)

# 3. 读取并处理 PCA 坐标数据 (.eigenvec)
pca_data <- fread(eigen_vec_file, header = FALSE)

# 合并 ID (SampleID) 并重命名 PC 列
pca_data <- pca_data %>%
  unite("SampleID", V1, V2, sep = "_", remove = TRUE) %>%
  rename_with(~ paste0("PC", seq_along(.)), .cols = -1)

# 4. 计算解释率 (使用 .eigenval 文件)
# 读取特征值
eigenvalues <- scan(eigen_val_file)

# 计算百分比 (特征值 / 总和)
pve <- eigenvalues / sum(eigenvalues) * 100

# 生成坐标轴标签文本 (保留2位小数)
x_lab_text <- paste0("PC", x_pc, " (", round(pve[x_pc], 2), "%)")
y_lab_text <- paste0("PC", y_pc, " (", round(pve[y_pc], 2), "%)")

# 5. 读取群体信息并合并
pop_data <- read.table(pop_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(pop_data) <- c("SampleID", "Group", "Region")

# 合并数据
plot_data <- inner_join(pca_data, pop_data, by = "SampleID")
print(paste("合并后保留样本数:", nrow(plot_data)))

# === 定义颜色和形状 ===
my_colors <- c("#f58073", "#edcb40", "#b4d66b", "#81b4d6") 
my_shapes <- c(18, 15, 16, 17) 

# 6. 绘图
p1 <- ggplot(data = plot_data, aes(x = .data[[paste0("PC", x_pc)]], 
                                   y = .data[[paste0("PC", y_pc)]], 
                                   color = Region,  
                                   fill = Region,   # 椭圆填充色
                                   shape = Region)) +
  # 添加置信椭圆 
  # stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.2, show.legend = FALSE) +
  # 添加散点
  geom_point(size = 3, alpha = 0.8) + 
  # 样式设置
  scale_color_manual(values = my_colors) + 
  scale_fill_manual(values = my_colors) +  
  scale_shape_manual(values = my_shapes) + 
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  # 使用正确计算的标签
  labs(x = x_lab_text, 
       y = y_lab_text,
       title = paste0("PCA Plot (PC", x_pc, " vs PC", y_pc, ")"))

# 7. 保存结果
print(p1)

## ggsave(filename = paste(outpre, "pc", x_pc, y_pc, "pdf", sep = "."), plot = p1, width = 10, height = 8)
## ggsave(filename = paste(outpre, "pc", x_pc, y_pc, "png", sep = "."), plot = p1, width = 10, height = 8, dpi = 300)
