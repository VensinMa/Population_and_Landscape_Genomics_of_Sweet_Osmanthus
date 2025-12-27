# ==============================================================================
# Sweet Osmanthus Gradient Forest Analysis
# 模块：RGB 三合一可视化 (Local, Forward, Reverse Genetic Offset)
# ==============================================================================

# 0. 环境设置与参数定义
# ==============================================================================
rm(list = ls()) # 清空环境
library(ggplot2)
library(data.table)
library(dplyr)
library(scales)
library(patchwork) 

# --- 参数设置 ---
work_dir      <- "/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/gradientForest_2025/2.5m"
target_ssp    <- "ssp585"       # 排放情景 ("ssp245", "ssp585")
target_period <- "2081-2100"    # 时间段   ("2041-2060", "2061-2080", "2081-2100")
cutoff_val    <- 0.95           # 归一化截断值 (去除极端值影响)

# 设置工作目录
if(dir.exists(work_dir)) setwd(work_dir) else stop("工作目录不存在: ", work_dir)
message("当前工作目录: ", getwd())

# ==============================================================================
# 1. 数据读取与预处理
# ==============================================================================
message(">>> 1. 正在读取数据...")

# 定义输入文件路径
# 假设目录结构为: Local_Genetic_Offset/Local_Genetic_Offset_ssp585_2081-2100.csv
f_local <- file.path("Local_Genetic_Offset",   paste0("Local_Genetic_Offset_", target_ssp, "_", target_period, ".csv"))
f_fwd   <- file.path("Forward_Genetic_Offset", paste0("Forward_Genetic_Offset_", target_ssp, "_", target_period, ".csv"))
f_rev   <- file.path("Reverse_Genetic_Offset", paste0("Reverse_Genetic_Offset_", target_ssp, "_", target_period, ".csv"))

# 检查文件是否存在
files_check <- c(Local=f_local, Forward=f_fwd, Reverse=f_rev)
if(any(!file.exists(files_check))) {
  missing <- names(files_check)[!file.exists(files_check)]
  stop("以下文件未找到: ", paste(missing, collapse=", "))
}

# 读取数据 (fread 速度快)
df_local <- fread(f_local)
df_fwd   <- fread(f_fwd)
df_rev   <- fread(f_rev)

# 重命名坐标列 (前两列是 lon/lat 或 x/y)
names(df_local)[1:2] <- c("x", "y")
names(df_fwd)[1:2]   <- c("x", "y")
names(df_rev)[1:2]   <- c("x", "y")

# 提取并重命名 Offset 列，准备合并
# 确保这里选择了正确的 Offset 列名 (通常是第3列或名为 offset 的列)
df_m <- df_local %>% dplyr::select(x, y, local_val = offset) %>%
  inner_join(df_fwd %>% dplyr::select(x, y, forward_val = offset), by = c("x", "y")) %>%
  inner_join(df_rev %>% dplyr::select(x, y, reverse_val = offset), by = c("x", "y"))

message("数据合并完成，共 ", nrow(df_m), " 个网格点。")

# ==============================================================================
# 2. RGB 数据归一化 (Normalization)
# ==============================================================================
message(">>> 2. 正在计算 RGB 颜色通道...")

# 定义归一化函数：使用分位数截断，防止极少数极大值导致整体颜色偏暗
rescale_channel <- function(x, quantile_prob = 0.95) {
  limit <- quantile(x, probs = quantile_prob, na.rm = TRUE)
  x_clipped <- ifelse(x > limit, limit, x)
  # Min-Max Normalization to [0, 1]
  return((x_clipped - min(x_clipped, na.rm=TRUE)) / (limit - min(x_clipped, na.rm=TRUE)))
}

df_rgb <- df_m %>%
  mutate(
    R = rescale_channel(local_val,   cutoff_val), # Red   = Local
    G = rescale_channel(forward_val, cutoff_val), # Green = Forward
    B = rescale_channel(reverse_val, cutoff_val), # Blue  = Reverse
    color_code = rgb(R, G, B)
  )

# ==============================================================================
# 3. 绘制 RGB 空间分布图 (Spatial Map)
# ==============================================================================
message(">>> 3. 绘制空间分布图...")

p_map <- ggplot(df_rgb, aes(x = x, y = y)) +
  geom_tile(aes(fill = color_code)) +
  scale_fill_identity() +
  coord_fixed() +
  theme_bw() +
  labs(
    title = paste0("Genetic Offset RGB Map (", target_ssp, " ", target_period, ")"),
    subtitle = paste0("Cutoff: ", cutoff_val*100, "% | R=Local, G=Forward, B=Reverse"),
    x = "Longitude", y = "Latitude"
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size=16, face="bold"),
    plot.subtitle = element_text(size=12)
  ) +
  # 添加左下角文字说明
  annotate("text", x = min(df_rgb$x), y = min(df_rgb$y), 
           label = "Red: Local Risk\nGreen: Migration Need\nBlue: Novel Environment", 
           hjust = 0, vjust = 0, color = "white", size = 3.5, fontface = "bold")
p_map
# 保存地图
out_map_name <- paste0("RGB_Map_", target_ssp, "_", target_period, "_cutoff", cutoff_val, ".pdf")
ggsave(out_map_name, p_map, width = 10, height = 8)
message("地图已保存: ", out_map_name)

# ==============================================================================
# 4. 绘制 2D 关联散点图 (Scatter Plots)
# ==============================================================================
message(">>> 4. 绘制关联散点图...")

# 定义通用的散点绘图函数
plot_rgb_scatter <- function(data, x_col, y_col, x_lab, y_lab) {
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(aes(color = color_code), size = 0.5, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.6, linetype = "dashed") +
    scale_color_identity() +
    labs(x = x_lab, y = y_lab) +
    theme_bw() +
    theme(
      panel.grid = element_line(color = "grey92"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    )
}

# 生成三张子图
p1 <- plot_rgb_scatter(df_rgb, "local_val",   "forward_val", "Local Offset",   "Forward Offset")
p2 <- plot_rgb_scatter(df_rgb, "forward_val", "reverse_val", "Forward Offset", "Reverse Offset")
p3 <- plot_rgb_scatter(df_rgb, "local_val",   "reverse_val", "Local Offset",   "Reverse Offset")

# 使用 patchwork 拼接
combined_scatter <- p1 | p2 | p3
combined_scatter
# 保存散点图
out_scatter_name <- paste0("RGB_Scatter_", target_ssp, "_", target_period, "_cutoff", cutoff_val, ".pdf")
ggsave(out_scatter_name, combined_scatter, width = 12, height = 4)
message("散点图已保存: ", out_scatter_name)

message("\n=== 所有任务完成 ===")