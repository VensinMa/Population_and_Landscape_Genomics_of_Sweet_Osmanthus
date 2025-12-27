# ==============================================================================
# 0. 环境设置与包加载
# ==============================================================================
# rm(list = ls()) # 清空环境变量
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
message("当前工作目录: ", getwd())

# 设置输入输出路径
input_climate_dir <- "extracted_data" 
Forward_Genetic_Offset_dir <- "Forward_Genetic_Offset" 
if (!dir.exists(Forward_Genetic_Offset_dir)) dir.create(Forward_Genetic_Offset_dir, recursive = TRUE)

# 定义未来时期列表
periods <- c("ssp245_2041-2060", "ssp245_2061-2080", "ssp245_2081-2100",
             "ssp585_2041-2060", "ssp585_2061-2080", "ssp585_2081-2100")

# 设置并行核心数 (根据你的服务器配置调整)
num_cores <- 28 

# ==============================================================================
# 1. 加载模型与确定变量
# ==============================================================================
model_file <- "gf.mod.10290.RData"
if(file.exists(model_file)){
  message("正在加载梯度森林模型: ", model_file)
  load(model_file)
} else {
  stop("错误: 找不到模型文件 '", model_file, "'")
}

if(!exists("gf.mod")) stop("模型文件加载后未找到 'gf.mod' 对象")
PredictEnvs <- names(gf.mod$X) 
message("模型使用的预测变量: ", paste(PredictEnvs, collapse = ", "))

# ==============================================================================
# 2. 处理当前气候数据 (Baseline)
# ==============================================================================
message("\n正在处理当前气候数据...")
current_file <- file.path(input_climate_dir, "Climate_current_66239_grids.csv")
if(!file.exists(current_file)) stop("找不到当前气候文件")

All_current_xy_envs <- fread(current_file)
All_current_xy_envs <- as.data.frame(All_current_xy_envs)
cols_needed <- c("lon", "lat", PredictEnvs)

# 筛选完整数据
All_current_xy_envs <- All_current_xy_envs[complete.cases(All_current_xy_envs[, cols_needed]), cols_needed]

# 转换当前数据到梯度森林生物空间 (Transformed)
current_pred_trans <- predict(gf.mod, All_current_xy_envs[, PredictEnvs])
current_coords <- All_current_xy_envs[, c("lon", "lat")]

message("当前数据处理完毕，有效样点数: ", nrow(current_pred_trans))

# ==============================================================================
# 3. 循环计算 Forward Genetic Offset
# ==============================================================================
message("\n开始计算 Forward Offset (严格匹配模式)...")

# 注册并行后端
cl <- makeCluster(num_cores)
registerDoParallel(cl)
message("已启动并行计算，核心数: ", num_cores)

for (period in periods) {
  message(">> 正在处理时期: ", period, " ...")
  
  # 3.1 读取未来数据
  file_name <- file.path(input_climate_dir, paste0("Climate_", period, "_66239_grids.csv"))
  if (!file.exists(file_name)) { warning("文件跳过: ", file_name); next }
  
  future_data <- fread(file_name)
  future_data <- as.data.frame(future_data)
  future_data <- future_data[complete.cases(future_data[, PredictEnvs]), cols_needed]
  
  # 3.2 转换未来数据到生物空间
  future_pred_trans <- predict(gf.mod, future_data[, PredictEnvs])
  future_coords <- future_data[, c("lon", "lat")]
  
  # 3.3 并行计算核心部分
  results_list <- foreach(i = 1:nrow(current_pred_trans), 
                          .packages = c("fields", "geosphere"), 
                          .combine = rbind) %dopar% {
                            
                            # A. 获取当前焦点种群的数据
                            focal_trans <- current_pred_trans[i, , drop=FALSE]
                            focal_coord <- current_coords[i, ]
                            
                            # B. 计算 Local Offset (原地偏移)
                            local_dist <- sqrt(sum((focal_trans - future_pred_trans[i, ])^2))
                            
                            # C. 计算 Forward Offset (全境搜索)
                            # fields::rdist 极其高效
                            bio_dists <- fields::rdist(focal_trans, future_pred_trans) 
                            
                            # 找到最小的生物距离
                            min_bio_dist <- min(bio_dists)
                            
                            # D. 筛选最佳候选点 (严格匹配，无容差)
                            candidates_idx <- which(bio_dists == min_bio_dist)
                            
                            # 获取候选点坐标
                            candidate_coords <- future_coords[candidates_idx, ]
                            
                            # 计算当前点到候选点的地理距离
                            geo_dists <- geosphere::distGeo(focal_coord, candidate_coords)
                            
                            # 在生物距离完全相等的候选点中，选择地理距离最近的
                            best_match_idx_in_candidates <- which.min(geo_dists)
                            best_match_global_idx <- candidates_idx[best_match_idx_in_candidates]
                            
                            # 获取最终指标
                            final_forward_offset <- min_bio_dist 
                            final_pred_dist <- geo_dists[best_match_idx_in_candidates] 
                            final_pred_dist_km <- final_pred_dist / 1000 
                            
                            # 计算方位角
                            best_coord <- future_coords[best_match_global_idx, ]
                            final_bearing <- geosphere::bearing(focal_coord, best_coord)
                            
                            # E. 返回结果
                            return(c(
                              x1 = as.numeric(focal_coord[1]),
                              y1 = as.numeric(focal_coord[2]),
                              local_offset = as.numeric(local_dist),
                              offset = as.numeric(final_forward_offset),
                              pred_dist_km = as.numeric(final_pred_dist_km),
                              bearing = as.numeric(final_bearing),
                              x2 = as.numeric(best_coord[1]),
                              y2 = as.numeric(best_coord[2])
                            ))
                          }
  
  # 3.4 转换结果并保存
  results_df <- as.data.frame(results_list)
  output_file <- file.path(Forward_Genetic_Offset_dir, paste0("Forward_Genetic_Offset_", period, ".csv"))
  fwrite(results_df, file = output_file, row.names = FALSE)
  
  message("  已保存: ", output_file)
  gc() # 内存清理
}

# 停止并行集群
stopCluster(cl)
message("\n所有任务完成！")

# ==============================================================================
# 4. 简单的可视化检查
# ==============================================================================
library(ggplot2)
library(viridis)

# 读取其中一个结果文件 (示例：ssp245_2081-2100)
test_period <- "ssp245_2081-2100"
test_file <- file.path(Forward_Genetic_Offset_dir, paste0("Forward_Genetic_Offset_", test_period, ".csv"))

if(file.exists(test_file)){
  plot_data <- read.csv(test_file)
  
  p <- ggplot(plot_data, aes(x = x1, y = y1, fill = offset)) +
    geom_tile() + # 如果是栅格点数据，tile 效果最好
    scale_fill_viridis(option = "magma", direction = -1) + # 使用 magma 色带，颜色越深代表偏移越大（风险越高）
    coord_fixed() +
    theme_minimal() +
    labs(title = paste("Forward Genetic Offset:", test_period),
         fill = "Offset",
         x = "Longitude", y = "Latitude")
  
  print(p)
  message("已绘制 ", test_period, " 的预览图")
}
