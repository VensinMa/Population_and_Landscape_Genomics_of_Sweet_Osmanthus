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
Reverse_Genetic_Offset_dir <- "Reverse_Genetic_Offset" # 输出目录改为 Reverse
if (!dir.exists(Reverse_Genetic_Offset_dir)) dir.create(Reverse_Genetic_Offset_dir, recursive = TRUE)

# 定义未来时期列表
periods <- c("ssp245_2041-2060", "ssp245_2061-2080", "ssp245_2081-2100",
             "ssp585_2041-2060", "ssp585_2061-2080", "ssp585_2081-2100")

# 设置并行核心数 (根据你的服务器配置调整，建议保留几颗核心给系统)
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
# 注意：这部分数据在后续循环中保持不变，作为“参考库”
current_pred_trans <- predict(gf.mod, All_current_xy_envs[, PredictEnvs])
current_coords <- All_current_xy_envs[, c("lon", "lat")]

message("当前数据处理完毕，有效样点数: ", nrow(current_pred_trans))

# ==============================================================================
# 3. 循环计算 Reverse Genetic Offset
# ==============================================================================
message("\n开始计算 Reverse Offset (这可能需要较长时间)...")

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
  
  # 检查行数一致性 (为了保证Local Offset也能同时算出来作为参考)
  if (nrow(future_data) != nrow(All_current_xy_envs)) {
    warning("行数不匹配，跳过: ", period); next 
  }
  
  # 3.2 转换未来数据到生物空间
  future_pred_trans <- predict(gf.mod, future_data[, PredictEnvs])
  future_coords <- future_data[, c("lon", "lat")]
  
  # 3.3 并行计算核心部分
  # 将大矩阵传递给 foreach 可能会有内存压力，但对于 6.6万行数据通常没问题
  # 我们需要遍历每一个“未来点” (i)，去“当前矩阵”里找最佳匹配
  
  results_list <- foreach(i = 1:nrow(future_pred_trans), 
                          .packages = c("fields", "geosphere"), 
                          .combine = rbind) %dopar% {
                            
                            # A. 获取当前焦点（未来时刻）的数据
                            focal_future_trans <- future_pred_trans[i, , drop=FALSE] # 必须保持矩阵格式
                            focal_future_coord <- future_coords[i, ]
                            
                            # B. 计算生物距离 (Bioclimatic Distance)
                            # 计算该未来点 i 到 所有当前点 的距离
                            # fields::rdist 效率非常高
                            dists <- fields::rdist(focal_future_trans, current_pred_trans) 
                            
                            # C. 提取 Local Offset (原地偏移)
                            # 即 dists 向量中第 i 个位置的值（对应地理位置相同的点）
                            local_val <- dists[i]
                            
                            # D. 提取 Reverse Offset (反向/最佳匹配偏移)
                            # 找到所有当前点中，生物环境与该未来点最相似的最小值
                            min_val <- min(dists)
                            
                            # E. 处理 Tie-breaking (多个最小值的情况)
                            # 找到所有具有最小值的索引
                            candidates_idx <- which(dists == min_val)
                            
                            # 如果有多个最佳匹配点，选地理距离最近的
                            if (length(candidates_idx) > 1) {
                              candidate_coords <- current_coords[candidates_idx, ]
                              geo_dists <- geosphere::distGeo(focal_future_coord, candidate_coords)
                              best_match_idx_local <- which.min(geo_dists) # 在候选者中的索引
                              best_match_idx_global <- candidates_idx[best_match_idx_local] # 全局索引
                              
                              final_pred_dist <- geo_dists[best_match_idx_local]
                            } else {
                              best_match_idx_global <- candidates_idx[1]
                              final_pred_dist <- geosphere::distGeo(focal_future_coord, current_coords[best_match_idx_global, ])
                            }
                            
                            # 获取最佳匹配点的坐标
                            best_coord <- current_coords[best_match_idx_global, ]
                            
                            # 计算方位角
                            final_bearing <- geosphere::bearing(focal_future_coord, best_coord)
                            
                            # F. 返回结果行
                            # 格式: x1, y1 (未来点坐标), Local_Offset, Reverse_Offset, Pred_Dist, Bearing, x2, y2 (当前最佳匹配坐标)
                            return(c(
                              x1 = as.numeric(focal_future_coord[1]),
                              y1 = as.numeric(focal_future_coord[2]),
                              local_offset = as.numeric(local_val),
                              offset = as.numeric(min_val),
                              pred_dist = as.numeric(final_pred_dist),
                              bearing = as.numeric(final_bearing),
                              x2 = as.numeric(best_coord[1]),
                              y2 = as.numeric(best_coord[2])
                            ))
                          }
  
  # 3.4 转换结果并保存
  results_df <- as.data.frame(results_list)
  output_file <- file.path(Reverse_Genetic_Offset_dir, paste0("Reverse_Genetic_Offset_", period, ".csv"))
  fwrite(results_df, file = output_file, row.names = FALSE)
  
  message("  已保存: ", output_file)
  # 显式清理内存
  gc()
}

# 停止并行集群
stopCluster(cl)
message("\n所有任务完成！")

# ==============================================================================
# 4. 可视化检查 Reverse Offset
# ==============================================================================
library(ggplot2)
library(viridis)

# 读取其中一个结果文件
test_period <- "ssp245_2081-2100"
test_file <- file.path(Reverse_Genetic_Offset_dir, paste0("Reverse_Genetic_Offset_", test_period, ".csv"))

if(file.exists(test_file)){
  plot_data <- fread(test_file)
  
  p <- ggplot(plot_data, aes(x = x1, y = y1, fill = offset)) +
    geom_tile() +
    scale_fill_viridis(option = "plasma", direction = -1) +
    coord_fixed() +
    theme_minimal() +
    labs(title = paste("Reverse Genetic Offset:", test_period),
         subtitle = "Minimum genetic distance to current analogue",
         fill = "Rev Offset")
  
  print(p)
  message("已绘制 ", test_period, " 的预览图")
}
