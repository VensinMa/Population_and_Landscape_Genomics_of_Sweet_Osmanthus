# ==============================================================================
# 0. 环境设置与包加载
# ==============================================================================
# rm(list = ls()) 
library(fields)
library(gdm)
library(geosphere)
library(gradientForest)
library(data.table)
library(doParallel)
library(foreach)

# 设置工作目录
workspace_dir <- "/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/gradientForest_2025/2.5m"
if (!dir.exists(workspace_dir)) dir.create(workspace_dir, recursive = TRUE)
setwd(workspace_dir)

# 输入输出路径
input_climate_dir <- "extracted_data" 
base_output_dir <- "Forward_Genetic_Offset_Limited" 
if (!dir.exists(base_output_dir)) dir.create(base_output_dir, recursive = TRUE)

# 定义时期
periods <- c("ssp245_2041-2060", "ssp245_2061-2080", "ssp245_2081-2100",
             "ssp585_2041-2060", "ssp585_2061-2080", "ssp585_2081-2100")

num_cores <- 28 

# ==============================================================================
# 1. 定义距离阈值
# ==============================================================================
dist_limits_km <- c(10, 20, 50, 100, 200, 500, 1000)
dist_limits_m  <- dist_limits_km * 1000

message("将分别生成以下距离限制的文件 (km): ", paste(dist_limits_km, collapse = ", "))

# ==============================================================================
# 2. 加载模型与 Baseline 数据
# ==============================================================================
model_file <- "gf.mod.10290.RData"
load(model_file)
PredictEnvs <- names(gf.mod$X) 

message("\n正在读取当前气候数据...")
current_file <- file.path(input_climate_dir, "Climate_current_66239_grids.csv")
All_current_xy_envs <- fread(current_file)
All_current_xy_envs <- as.data.frame(All_current_xy_envs)
cols_needed <- c("lon", "lat", PredictEnvs)
All_current_xy_envs <- All_current_xy_envs[complete.cases(All_current_xy_envs[, cols_needed]), cols_needed]

# 转换当前数据
current_pred_trans <- predict(gf.mod, All_current_xy_envs[, PredictEnvs])
current_coords <- All_current_xy_envs[, c("lon", "lat")]

message("当前数据准备完毕。")

# ==============================================================================
# 3. 循环计算 (并行核心部分)
# ==============================================================================
cl <- makeCluster(num_cores)
registerDoParallel(cl)

for (period in periods) {
  message("\n>> 正在处理时期: ", period)
  
  # 3.1 读取未来数据
  file_name <- file.path(input_climate_dir, paste0("Climate_", period, "_66239_grids.csv"))
  if (!file.exists(file_name)) next
  
  future_data <- fread(file_name)
  future_data <- as.data.frame(future_data)
  future_data <- future_data[complete.cases(future_data[, PredictEnvs]), cols_needed]
  
  if (nrow(future_data) != nrow(All_current_xy_envs)) { warning("行数不匹配，跳过"); next }
  
  # 3.2 转换未来数据
  future_pred_trans <- predict(gf.mod, future_data[, PredictEnvs])
  future_coords <- future_data[, c("lon", "lat")]
  
  # 3.3 并行计算
  message("   正在并行计算各类距离限制下的 Offset...")
  
  results_matrix <- foreach(i = 1:nrow(current_pred_trans), 
                            .packages = c("fields", "geosphere"), 
                            .combine = rbind) %dopar% {
                              
                              # [数据准备]
                              focal_trans <- current_pred_trans[i, , drop=FALSE]
                              focal_coord <- current_coords[i, ]
                              
                              # 计算 Local Offset
                              local_offset <- sqrt(sum((focal_trans - future_pred_trans[i, ])^2))
                              
                              # 预计算：该点到未来所有点的 生物距离 和 地理距离
                              all_gfOffset <- as.vector(fields::rdist(focal_trans, future_pred_trans))
                              all_geoDists <- geosphere::distGeo(focal_coord, future_coords)
                              
                              # 初始化输出向量：前3列是基础信息
                              row_out <- c(as.numeric(focal_coord[1]), as.numeric(focal_coord[2]), as.numeric(local_offset))
                              
                              # [核心逻辑] 遍历每一个距离限制
                              for (limit_m in dist_limits_m) {
                                
                                # 1. 距离筛选 (Filtering)
                                valid_indices <- which(all_geoDists < limit_m)
                                # 这里移除了 length(valid_indices)==0 的判断，默认一定会找到点
                                
                                # 2. 提取有效数据
                                valid_gfOffset <- all_gfOffset[valid_indices]
                                valid_geoDists <- all_geoDists[valid_indices]
                                valid_indices_global <- valid_indices 
                                
                                # 3. 找最小生物距离
                                min_gf_val <- min(valid_gfOffset)
                                best_bio_indices <- which(valid_gfOffset == min_gf_val)
                                
                                # 4. 如果并列，找最小地理距离 (Tie-breaking)
                                candidates_geo <- valid_geoDists[best_bio_indices]
                                min_geo_val <- min(candidates_geo)
                                best_geo_idx_local <- which(candidates_geo == min_geo_val)
                                
                                # 5. 如果还有并列，随机抽样
                                if (length(best_geo_idx_local) > 1) {
                                  final_pick_idx <- sample(best_geo_idx_local, 1)
                                } else {
                                  final_pick_idx <- best_geo_idx_local[1]
                                }
                                
                                # 6. 追溯坐标
                                final_valid_idx_local <- best_bio_indices[final_pick_idx]
                                best_global_idx <- valid_indices_global[final_valid_idx_local]
                                
                                best_pt_coord <- future_coords[best_global_idx, ]
                                best_bearing <- geosphere::bearing(focal_coord, best_pt_coord)
                                
                                # 追加结果
                                row_out <- c(row_out, 
                                             as.numeric(min_gf_val),      # forwardOffset
                                             as.numeric(min_geo_val),     # predDist (meters)
                                             as.numeric(best_bearing),    # bearing
                                             as.numeric(best_pt_coord[1]), # x2
                                             as.numeric(best_pt_coord[2])) # y2
                              }
                              return(row_out)
                            }
  
  # 3.4 拆分并保存
  message("   正在拆分并保存 CSV 文件...")
  results_df <- as.data.frame(results_matrix)
  
  n_metrics <- 5 
  
  for (k in 1:length(dist_limits_km)) {
    km_val <- dist_limits_km[k]
    start_col <- 3 + (k - 1) * n_metrics + 1
    end_col   <- 3 + k * n_metrics
    
    sub_df <- results_df[, c(1:3, start_col:end_col)]
    colnames(sub_df) <- c("x1", "y1", "local_offset", 
                          "offset", "pred_dist_m", "bearing", "x2", "y2")
    
    file_out_name <- file.path(base_output_dir, 
                               paste0("Forward_Genetic_Offset_", km_val, "km_", period, ".csv"))
    
    fwrite(sub_df, file_out_name, row.names = FALSE)
  }
  
  message("   时期 ", period, " 完成。")
  gc() 
}

stopCluster(cl)
message("\n所有任务完成！")



# ==============================================================================
# 4. 绘图
# ==============================================================================
# ==============================================================================
# 1. 环境与数据加载
# ==============================================================================
library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)

# 设置为你当前的 Base 工作目录 (保持不变)
# setwd("/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/gradientForest_2025/2.5m")

# 定义存放 CSV 的子文件夹名称
input_subdir <- "Forward_Genetic_Offset_Limited"

# [修改点1] 定义所有三个时间段
timeranges <- c("2041-2060", "2061-2080", "2081-2100")

# 定义文件名前缀
file_prefixes <- c(
  "Forward_Genetic_Offset_10km_ssp245_",   "Forward_Genetic_Offset_10km_ssp585_",
  "Forward_Genetic_Offset_20km_ssp245_",   "Forward_Genetic_Offset_20km_ssp585_",
  "Forward_Genetic_Offset_50km_ssp245_",   "Forward_Genetic_Offset_50km_ssp585_",
  "Forward_Genetic_Offset_100km_ssp245_",  "Forward_Genetic_Offset_100km_ssp585_",
  "Forward_Genetic_Offset_200km_ssp245_",  "Forward_Genetic_Offset_200km_ssp585_",
  "Forward_Genetic_Offset_500km_ssp245_",  "Forward_Genetic_Offset_500km_ssp585_",
  "Forward_Genetic_Offset_1000km_ssp245_", "Forward_Genetic_Offset_1000km_ssp585_"
)

# [修改点2] 生成所有时间段的所有文件路径组合
# 使用 expand.grid 创建所有组合，然后拼接
combinations <- expand.grid(prefix = file_prefixes, time = timeranges)
file_names <- paste0(combinations$prefix, combinations$time, ".csv")
file_paths <- file.path(input_subdir, file_names)

# 初始化数据列表
data_list <- list()

# 遍历读取
for (file in file_paths) {
  if (!file.exists(file)) {
    # 仅在文件缺失时警告，不中断
    # warning(paste("文件不存在，跳过:", file))
    next
  }
  
  # 读取数据
  df <- fread(file) 
  
  # --- 文件名解析 ---
  base_name <- basename(file)
  file_info <- strsplit(base_name, "_")[[1]]
  
  # 提取距离 (第4个元素)
  dist_str <- file_info[4] 
  dist_val <- as.numeric(gsub("km", "", dist_str))
  
  # 提取情景 (第5个元素)
  scenario_str <- file_info[5]
  
  # 提取时间 (第6个元素)
  time_str <- gsub(".csv", "", file_info[6])
  
  # 添加新列
  df$search_distance <- dist_val
  df$scenario <- scenario_str
  df$time_range <- time_str
  
  data_list[[file]] <- df
}

# 合并所有数据
all_data <- rbindlist(data_list, use.names = TRUE, fill = TRUE)

# [修改点3] 统一 Offset 列名 (兼容 forward_offset, offset, forwardOffset)
if ("forward_offset" %in% names(all_data)) {
  all_data$forwardOffset <- all_data$forward_offset
} else if ("offset" %in% names(all_data)) {
  all_data$forwardOffset <- all_data$offset
}

# [修改点4] 设置时间因子的顺序，确保画图时按时间排序
all_data$time_range <- factor(all_data$time_range, levels = c("2041-2060", "2061-2080", "2081-2100"))

# ==============================================================================
# 2. 数据汇总
# ==============================================================================
summary_data <- all_data %>%
  group_by(search_distance, scenario, time_range) %>%
  summarise(
    median_offset = median(forwardOffset, na.rm = TRUE), # 确保使用统一后的列名
    lower_offset = quantile(forwardOffset, 0.25, na.rm = TRUE),
    upper_offset = quantile(forwardOffset, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# ==============================================================================
# 3. 绘图 (分面展示)
# ==============================================================================


p <- ggplot(summary_data, aes(x = search_distance, y = median_offset, 
                              group = scenario, color = scenario, fill = scenario)) +
  # 绘制置信区间
  geom_ribbon(aes(ymin = lower_offset, ymax = upper_offset), alpha = 0.2, color = NA) +
  # 绘制折线
  geom_line(linewidth = 0.8) +
  # 绘制点
  geom_point(size = 2, shape = 21, color = "white", stroke = 0.8) + 
  
  # [修改点5] 使用 facet_wrap 分面显示三个时间段
  facet_wrap(~ time_range, ncol = 3) + 
  
  # 坐标轴与标签
  labs(
    title = "Forward Genetic Offset via Dispersal Distance",
    subtitle = "Comparison across SSP scenarios and time periods",
    x = "Dispersal Distance Limit (km)",
    y = "Forward Genetic Offset",
    color = "Scenario",
    fill = "Scenario"
  ) +
  
  # 颜色设置
  scale_color_manual(values = c("ssp245" = "#4B5CC4", "ssp585" = "#D6404E")) +
  scale_fill_manual(values = c("ssp245" = "#4B5CC4", "ssp585" = "#D6404E")) +
  
  # X轴刻度
  scale_x_continuous(breaks = unique(summary_data$search_distance)) +
  
  # 主题
  theme_bw() +
  theme(
    # 图例位置调整到下方，因为横向三个图比较宽
    legend.position = "bottom",
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1), # X轴文字稍微倾斜防止重叠
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    
    # 分面标题样式
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 11, face = "bold"),
    
    panel.grid.minor = element_blank()
  )

print(p)

# 保存 (图片稍微宽一点，以容纳三个分面)
ggsave("Forward_Offset_Curve_All_Periods.pdf", p, width = 12, height = 5)
message("图片已保存至当前目录: ", getwd(), "/Forward_Offset_Curve_All_Periods.pdf")