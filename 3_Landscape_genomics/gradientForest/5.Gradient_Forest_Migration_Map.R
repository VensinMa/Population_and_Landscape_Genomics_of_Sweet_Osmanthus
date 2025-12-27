# ==============================================================================
# 绘制 Forward Genetic Offset 的迁移距离、方向及极坐标直方图
# (使用自定义 SHP 底图版)
# ==============================================================================

rm(list = ls()) # 清空环境
library(ggplot2)
library(data.table)
library(sf)            # 处理 SHP 文件
library(ggspatial)     # 用于添加指北针和比例尺
library(scales)
library(patchwork)     # 用于拼图

# ==============================================================================
# 1. 参数设置
# ==============================================================================
work_dir      <- "/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/gradientForest_2025/2.5m"
shp_path      <- "/home/data/地图/China map data/国界_4m/bou1_4p.shp" # SHP 文件路径
target_ssp    <- "ssp585"       # 排放情景
target_period <- "2081-2100"    # 时间段
cutoff_val    <- 0.95           # 用于确定距离显示的上限

# 设置工作目录
if(dir.exists(work_dir)) setwd(work_dir) else stop("工作目录不存在: ", work_dir)

# 构建输入文件路径
input_dir <- "Forward_Genetic_Offset"
file_name <- paste0("Forward_Genetic_Offset_", target_ssp, "_", target_period, ".csv")
file_path <- file.path(input_dir, file_name)

# ==============================================================================
# 2. 读取与预处理数据
# ==============================================================================
# 2.1 读取 CSV 数据
message("正在读取数据文件: ", file_path)
if(!file.exists(file_path)) stop("CSV 文件不存在: ", file_path)
df <- fread(file_path)

# 检查必要列
required_cols <- c("x1", "y1", "pred_dist_km", "bearing")
if(!all(required_cols %in% colnames(df))) stop("数据缺失必要的列!")

# --- 方位角转换 [-180, 180] -> [0, 360] ---
if (min(df$bearing, na.rm = TRUE) < 0) {
  message("检测到负值方位角，正在转换为 [0, 360] 范围...")
  df$bearing <- ifelse(df$bearing < 0, df$bearing + 360, df$bearing)
}

# 2.2 读取 SHP 底图数据
message("正在读取地图 SHP 文件: ", shp_path)
if(!file.exists(shp_path)) stop("SHP 文件不存在: ", shp_path)

# 读取 SHP
china_shp <- st_read(shp_path, quiet = TRUE)

# 检查并处理坐标系
current_crs <- st_crs(china_shp)

if (is.na(current_crs)) {
  # 情况 A: SHP 文件完全缺失坐标系信息 (常见于 bou1_4p)
  message("警告: SHP 文件缺失坐标系定义 (.prj)。假设其为 WGS84 并手动赋值...")
  st_crs(china_shp) <- 4326
} else {
  # 情况 B: 有坐标系，强制转换为 WGS84 以匹配气象数据
  # (如果原本就是 4326，这一步是安全的，不会改变数据)
  message("正在统一底图坐标系为 WGS84...")
  china_shp <- st_transform(china_shp, 4326)
}
# ==============================================================================
# 3. 准备极坐标直方图 (Polar Histogram / Wind Rose)
# ==============================================================================
message("正在生成极坐标直方图...")

bin_width <- 3
breaks_seq <- seq(0, 360, by = bin_width)
bearing_cut <- cut(df$bearing, breaks = breaks_seq, include.lowest = TRUE, right = FALSE)
bearing_counts <- as.data.frame(table(bearing_cut))
bearing_counts$mid_angle <- head(breaks_seq, -1) + bin_width / 2
bearing_counts$count <- bearing_counts$Freq

cyclic_colors <- c(
  "#698e6a", # N (0)
  "#69559b", # E (90)  
  "#af5944", # S (180) 
  "#f6fc71", # W (270) 
  "#698e6a"  # N (360) 
)

p_rose <- ggplot(bearing_counts, aes(x = mid_angle, y = count, fill = mid_angle)) +
  geom_bar(stat = "identity", width = bin_width, color = "grey30", linewidth = 0) +
  scale_x_continuous(limits = c(0, 360), breaks = c(0, 90, 180, 270), labels = c("N", "E", "S", "W")) +
  scale_fill_gradientn(colors = cyclic_colors, guide = "none") +
  coord_polar(start = 0) +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(
    #axis.text.y = element_blank(), 不显示方位和数量的标注
    axis.text.y = element_text(size = 5, face = "bold", color = "black"),
    axis.text.x = element_text(size = 8, face = "bold", color = "black"),
    panel.grid.major = element_line(color = "grey80", size = 0.3),
    plot.background = element_rect(fill = "transparent", color = NA)
  )
p_rose # 现在看着有点小，但后边叠加到主图上就好了

# ==============================================================================
# 4. 绘图 (e): 迁移距离 (Distance km)
# ==============================================================================
message("正在绘制迁移距离图...")
limit_dist <- quantile(df$pred_dist_km, probs = cutoff_val, na.rm = TRUE)
df$plot_dist <- ifelse(df$pred_dist_km > limit_dist, limit_dist, df$pred_dist_km)
# ==============================================================================
# 最终绘图代码 (a) Migration Distance
# ==============================================================================
p_dist <- ggplot() +
  # 使用 borders 绘制底图
  # 注意！borders 函数绘制的中国地图缺少台湾省，建议使用自己本地带审图号的shp
  # borders("world", regions = "China", fill = NA, colour = "grey50", linewidth = 0.2) +
  
  # 1. 添加 SHP 底图 (使用自定义填充颜色)
  geom_sf(data = china_shp, fill = "#e1e1e1", color = "#e1e1e1", size = 0.2) +
  
  # 2. 绘制栅格数据
  geom_tile(data = df, aes(x = x1, y = y1, fill = plot_dist)) +
  scale_fill_gradientn(
    colors = c("#3288bd", "#66c2a5", "#abdda4", "#e6f598", "#fee08b", "#fdae61", "#f46d43", "#d53e4f"),
    name = "Distance (km)",
    limits = c(0, limit_dist),
    oob = scales::squish
  ) +
  
  # 3. 设置坐标系 (自动裁剪到数据范围)
  coord_sf(
    xlim = c(min(df$x1), max(df$x1)), 
    ylim = c(min(df$y1), max(df$y1)), 
    crs = 4326
  ) +
  
  # 4. 指北针和比例尺
  annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  
  # 5. 主题与字体设置
  theme_bw() +
  labs(title = paste0("(a) Migration Distance: ", target_ssp, " ", target_period), 
       x = "Longitude", y = "Latitude") +
  theme(
    panel.grid = element_blank(),
    
    # 坐标轴刻度字体 (黑, 12号, 加粗)
    axis.text = element_text(color = "black", size = 12, face = "bold"),
    
    # 坐标轴标题字体 (黑, 14号, 加粗)
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    
    # 图例位置与背景
    legend.position = c(0.08, 0.82),
    legend.background = element_rect(fill = "white", color = "black", size = 0.1),
    
    # 图例文字大小 (使图例整体看起来更小更精致)
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.margin = margin(3, 3, 3, 3)
  ) +
  
  # 6. 精细控制图例色条尺寸 (使其变细变短)
  guides(fill = guide_colorbar(
    barwidth = 0.8,       # 宽度变细
    barheight = 4,        # 高度变短
    ticks.linewidth = 0.5,# 刻度线变细
    title.position = "top",
    title.hjust = 0.5
  ))

# 打印预览
print(p_dist)

# ==============================================================================
# 5. 绘图 (f): 迁移方向 (Bearing) + 嵌入风玫瑰图
# ==============================================================================
message("正在绘制迁移方向图...")

# 主地图
p_bear_map <- ggplot() +
  # 添加 SHP 底图
  geom_sf(data = china_shp, fill = "#e1e1e1", color = "#e1e1e1", size = 0.2) +
  
  # 注意borders函数绘制的中国地图缺少台湾省，建议使用带审图号的shp
  # borders("world", regions = "China", fill = NA, colour = "grey50", linewidth = 0.2) +
  
  geom_tile(data = df, aes(x = x1, y = y1, fill = bearing)) +
  
  scale_fill_gradientn(
    colors = cyclic_colors,
    name = "Bearing (°)",
    limits = c(0, 360),        # 因为会绘制极方图并叠加，所以这个图例不需要
    breaks = c(0, 90, 180, 270, 360),
    labels = c("N", "E", "S", "W", "N")
  ) +

  coord_sf(xlim = c(min(df$x1), max(df$x1)), ylim = c(min(df$y1), max(df$y1)), crs = 4326) +
  annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme_bw() +
  labs(title = paste0("(b) Migration Direction: ", target_ssp, " ", target_period), x = "Longitude", y = "Latitude") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 12, face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    legend.position = "none" # legend.position = "none" 关闭图例，"right"则在主图右侧显示
  )
p_bear_map

# 插入玫瑰图
p_bear_combined <- p_bear_map + 
  inset_element(p_rose, left = 0.01, bottom = 0.60, right = 0.35, top = 0.98, align_to = 'plot')
p_bear_combined
# ==============================================================================
# 6. 保存图片
# ==============================================================================
combined_plot <- p_dist | p_bear_combined
out_filename_base <- paste0("Migration_Map_With_Rose_", target_ssp, "_", target_period)

ggsave(paste0(out_filename_base, ".pdf"), combined_plot, width = 20, height = 5)
message("绘图完成！文件已保存: ", out_filename_base)
