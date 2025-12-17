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
eigen_file <- "202_samples_gcta_pca_result.eigenvec" # PCA 结果文件
pop_file   <- "202samples.pop"                       # 分组文件

# 输出设置
outpre <- "PCA_out_figure"
x_pc <- 1 # X轴使用第几主成分 (PC1)
y_pc <- 2 # Y轴使用第几主成分 (PC2)


# 3. 读取并处理 PCA 数据
# header = FALSE, GCTA输出通常没有表头
pca_data <- fread(eigen_file, header = FALSE)

# 合并前两列为 SampleID (对应你之前的需求)
# 注意：确保合并后的名字(如 CP_1)能和 pop 文件对应上
pca_data <- pca_data %>%
  unite("SampleID", V1, V2, sep = "_", remove = TRUE) %>%
  # 重命名列：第一列是ID，后面依次是 PC1, PC2...
  rename_with(~ paste0("PC", seq_along(.)), .cols = -1)

data = pca_data[2:4]

pca <- prcomp(data)  # data 为样本 × 变量矩阵
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
var_explained[1:3]

eigenval <- scan("202_samples_gcta_pca_result.eigenval")

var_explained <- eigenval / sum(eigenval)

round(var_explained[1:3] * 100, 2)



# 4. 读取群体信息文件
# 假设文件没有表头，第一列是样本名，第二列是群体(CP/CX)，第三列是区域(East)
pop_data <- read.table(pop_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# 给pop数据加表头方便后续操作
colnames(pop_data) <- c("SampleID", "Group", "Region")

# 5. 合并数据 (关键步骤)
# 使用 inner_join 自动根据 SampleID 对齐数据，防止顺序错乱
plot_data <- inner_join(pca_data, pop_data, by = "SampleID")

# 检查一下是否成功合并，如果行数变少说明ID没对上
print(paste("合并后保留样本数:", nrow(plot_data)))

# === 定义颜色和形状 ===
my_colors <- c("#f58073", "#edcb40", "#b4d66b", "#81b4d6") 
my_shapes <- c(18, 15, 16, 17) 

# === 计算方差用于标签 (保持不变) ===
pc_columns <- grep("^PC", colnames(plot_data))
x_lab_text <- paste0("PC", x_pc, " (", round(var(plot_data[[paste0("PC", x_pc)]]) / sum(apply(plot_data[, pc_columns], 2, var)) * 100, 2), "%)")
y_lab_text <- paste0("PC", y_pc, " (", round(var(plot_data[[paste0("PC", y_pc)]]) / sum(apply(plot_data[, pc_columns], 2, var)) * 100, 2), "%)")

# 6. 绘图
p1 <- ggplot(data = plot_data, aes(x = .data[[paste0("PC", x_pc)]], 
                                   y = .data[[paste0("PC", y_pc)]], 
                                   color = Region,  # 点的轮廓颜色
                                   fill = Region,   # 【新增】椭圆的填充颜色
                                   shape = Region)) +
  
  # --- 【新增】添加置信椭圆 ---
  # level = 0.95 表示 95% 置信区间
  # alpha = 0.2 设置透明度，防止挡住点
  # show.legend = FALSE 不在图例中显示椭圆的方块，保持图例整洁
  # stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.2, show.legend = FALSE) +
  # 添加点 (建议放在椭圆后面，这样点会浮在椭圆上方)
  geom_point(size = 3, alpha = 0.8) + 
  scale_color_manual(values = my_colors) + # 设置点的颜色
  scale_fill_manual(values = my_colors) +  # 【新增】设置椭圆填充颜色 (必须与color一致)
  scale_shape_manual(values = my_shapes) + # 设置形状
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  labs(x = x_lab_text, 
       y = y_lab_text,
       title = paste0("PCA Plot (PC", x_pc, " vs PC", y_pc, ")"))

# 7. 保存结果
print(p1)

#ggsave(filename = paste(outpre, "pc", x_pc, y_pc, "pdf", sep = "."), plot = p1, width = 10, height = 8)
#ggsave(filename = paste(outpre, "pc", x_pc, y_pc, "png", sep = "."), plot = p1, width = 10, height = 8, dpi = 300)
