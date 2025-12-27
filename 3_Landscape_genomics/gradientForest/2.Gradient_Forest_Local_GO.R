# ==============================================================================
# 0. 环境设置与包加载
# ==============================================================================
# rm(list = ls()) # 清空环境变量
library(fields)
library(gdm)
library(geosphere)
library(gradientForest)
library(data.table)

# 设置工作目录
workspace_dir <- "/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/gradientForest_2025/2.5m"
if (!dir.exists(workspace_dir)) dir.create(workspace_dir, recursive = TRUE)
setwd(workspace_dir)
message("当前工作目录: ", getwd())

# 设置输入输出路径
input_climate_dir <- "extracted_data"    # 存放气候数据的文件夹（包括当前及未来气候情景）
Local_Genetic_Offset_dir <- "Local_Genetic_Offset"    # 存放结果的文件夹
if (!dir.exists(Local_Genetic_Offset_dir)) dir.create(Local_Genetic_Offset_dir, recursive = TRUE)

# 定义未来时期列表
periods <- c("ssp245_2041-2060", "ssp245_2061-2080", "ssp245_2081-2100",
             "ssp585_2041-2060", "ssp585_2061-2080", "ssp585_2081-2100")

# ==============================================================================
# 1. 加载模型与确定变量
# ==============================================================================
model_file <- "gf.mod.10290.RData"

if(file.exists(model_file)){
  message("正在加载梯度森林模型: ", model_file)
  load(model_file)
} else {
  stop("错误: 找不到模型文件 '", model_file, "'，请检查路径！")
}

# 从模型中提取所需的预测变量名
# 模型对象名为 gf.mod
if(!exists("gf.mod")) stop("模型文件加载后未找到 'gf.mod' 对象")
PredictEnvs <- names(gf.mod$X) 
message("模型使用的预测变量: ", paste(PredictEnvs, collapse = ", "))

# ==============================================================================
# 2. 处理当前气候数据 (Baseline)
# ==============================================================================
message("\n正在处理当前气候数据...")

# 读取当前气候数据
current_file <- file.path(input_climate_dir, "Climate_current_66239_grids.csv")
if(!file.exists(current_file)) stop("找不到当前气候文件: ", current_file)

All_current_xy_envs <- fread(current_file)
All_current_xy_envs <- as.data.frame(All_current_xy_envs)

# 筛选完整数据（经纬度 + 预测变量）
cols_needed <- c("lon", "lat", PredictEnvs)
# 检查列是否存在
missing_cols <- setdiff(cols_needed, names(All_current_xy_envs))
if(length(missing_cols) > 0) stop("当前气候数据缺失以下列: ", paste(missing_cols, collapse=", "))

All_current_xy_envs <- All_current_xy_envs[complete.cases(All_current_xy_envs[, cols_needed]), cols_needed]

# 计算当前的梯度森林转换值 (Transformed Importance)
# 结果只保留转换后的生物空间值
current_pred_trans <- predict(gf.mod, All_current_xy_envs[, PredictEnvs])
dim_current <- dim(current_pred_trans)
message("当前数据处理完毕，有效样点数: ", dim_current[1])

# ==============================================================================
# 3. 循环计算未来遗传偏移 (Genetic Offset)
# ==============================================================================
message("\n开始计算未来情景的遗传偏移...")

for (period in periods) {
  message("正在处理时期: ", period, " ...")
  
  # 3.1 读取未来数据
  file_name <- file.path(input_climate_dir, paste0("Climate_", period, "_66239_grids.csv"))
  
  if (!file.exists(file_name)) {
    warning("文件不存在，跳过: ", file_name)
    next
  }
  
  future_data <- fread(file_name)
  future_data <- as.data.frame(future_data)
  
  # 3.2 筛选与匹配
  # 注意：计算 Offset 必须保证未来数据的行数、顺序与当前数据完全一致
  # 假设 csv 文件的行顺序是对应的。如果文件包含不同样点，需要先 merge/join。
  future_data <- future_data[complete.cases(future_data[, PredictEnvs]), cols_needed]
  
  # [关键检查] 检查行数是否匹配
  if (nrow(future_data) != nrow(All_current_xy_envs)) {
    warning(paste0("警告: ", period, " 的数据行数 (", nrow(future_data), 
                   ") 与当前数据 (", nrow(All_current_xy_envs), ") 不一致！无法计算一一对应的偏移。跳过此文件。"))
    next
  }
  
  # 3.3 预测未来环境的转换值
  future_pred_trans <- predict(gf.mod, future_data[, PredictEnvs])
  
  # 3.4 计算欧几里得距离 (Genetic Offset)
  # 公式: sqrt( sum( (Future_trans - Current_trans)^2 ) )
  # 因为 predict 返回的是矩阵/数据框，可以直接相减
  diff_sq <- (future_pred_trans - current_pred_trans) ^ 2
  offsetAll <- sqrt(rowSums(diff_sq))
  
  # 3.5 整合结果
  Offset <- cbind(future_data[, c("lon", "lat")], offset = offsetAll)
  
  # 3.6 保存结果
  Local_Genetic_Offset_file <- file.path(Local_Genetic_Offset_dir, paste0("Local_Genetic_Offset_", period, ".csv"))
  fwrite(Offset, file = Local_Genetic_Offset_file, row.names = FALSE)
  
  message("  已保存: ", Local_Genetic_Offset_file)
}

message("\n所有任务完成！")

# ==============================================================================
# 4. 简单的可视化检查 (Visual Check)
# ==============================================================================
library(ggplot2)
library(viridis)

# 读取其中一个结果文件 (示例：ssp245_2081-2100)
test_period <- "ssp245_2081-2100"
test_file <- file.path(Local_Genetic_Offset_dir, paste0("Local_Genetic_Offset_", test_period, ".csv"))

if(file.exists(test_file)){
  plot_data <- read.csv(test_file)
  
  p <- ggplot(plot_data, aes(x = lon, y = lat, fill = offset)) +
    geom_tile() + # 如果是栅格点数据，tile 效果最好
    scale_fill_viridis(option = "magma", direction = -1) + # 使用 magma 色带，颜色越深代表偏移越大（风险越高）
    coord_fixed() +
    theme_minimal() +
    labs(title = paste("Local Genetic Offset:", test_period),
         fill = "Offset",
         x = "Longitude", y = "Latitude")
  
  print(p)
  message("已绘制 ", test_period, " 的预览图")
}
