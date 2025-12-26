library(terra)

# 设置工作目录
setwd("/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/gradientForest_2025/2.5m/")

# 读取坐标文件
coordinates <- read.csv("input/OF_AREA_2.5M.csv")
# 确保经纬度列为数值型
coordinates$lon <- as.numeric(coordinates$lon)
coordinates$lat <- as.numeric(coordinates$lat)
coords_mat <- coordinates[, c("lon", "lat")] # 提取坐标矩阵用于 terra

# 环境图层目录
tif_directory <- "/home/data/MIROC6_30s/future_climate"

# 定义情景和时间范围
scenario_time_combinations <- c("current",
                                "ssp245_2041-2060", "ssp245_2061-2080", "ssp245_2081-2100",
                                "ssp585_2041-2060", "ssp585_2061-2080", "ssp585_2081-2100"
)

# 初始化基础数据框 (仅保留样品名、经纬度，不再提取固定变量)
base_data <- coordinates[, c(1, 2, 3)] # 假设前三列是 ind, lon, lat
colnames(base_data)[1:3] <- c("ind", "lon", "lat")

# 2. 循环处理不同情景的气候变量 ----------------------------------------------------
output_directory <- "extracted_data"
if (!dir.exists(output_directory)) dir.create(output_directory)

for (combo in scenario_time_combinations) {
  cat(sprintf("正在处理情景: %s ...\n", combo))
  
  # 构建 BIO1-BIO19 的文件路径列表
  bio_files <- c()
  
  # 注意：Current 数据的文件名格式通常与未来情景不同，请根据实际情况修改下面的 pattern
  if (combo == "current") {
    # 假设 Current 文件名为: wc2.1_30s_bio_1.tif
    # 如果你的 current 文件名也带有 MIROC6，请修改此处
    for (i in 1:19) {
      # 示例 pattern: wc2.1_30s_bio_1.tif
      fname <- sprintf("wc2.1_30s_bio_%d.tif", i) 
      # 如果你的文件名是 wc2.1_30s_bioc_MIROC6_current_BIO1.tif，请取消下面这行的注释并修改
      fname <- sprintf("wc2.1_30s_bioc_MIROC6_current_BIO%d.tif", i)
      bio_files <- c(bio_files, file.path(tif_directory, fname))
    }
  } else {
    # 未来情景文件名
    for (i in 1:19) {
      # 示例 pattern: wc2.1_30s_bioc_MIROC6_ssp245_2041-2060_BIO1.tif
      fname <- sprintf("wc2.1_30s_bioc_MIROC6_%s_BIO%d.tif", combo, i)
      bio_files <- c(bio_files, file.path(tif_directory, fname))
    }
  }
  
  # 检查所有文件是否存在
  if (all(file.exists(bio_files))) {
    # 加载所有 BIO 图层为一个 Stack
    bio_stack <- rast(bio_files)
    names(bio_stack) <- paste0("BIO", 1:19) # 重命名图层以便识别
    
    # 一次性提取所有 19 个变量
    extracted_bios <- extract(bio_stack, coords_mat)
    
    # 合并数据 (extracted_bios 第一列是 ID，需去掉)
    final_result <- cbind(base_data, extracted_bios[, -1])
    
    # 3. 调整列顺序 ------------------------------------------------------------
    # 确保列名存在再排序，防止出错
    target_cols <- c("ind", "lon", "lat", paste0("BIO", 1:19))
    
    # 只选择存在的列
    cols_to_keep <- intersect(target_cols, colnames(final_result))
    final_result <- final_result[, cols_to_keep]
    
    # 4. 保存结果 --------------------------------------------------------------
    # 使用 202samples 以匹配输入文件，或者使用 nrow(final_result) 动态获取
    output_file_name <- sprintf("Climate_%s_%d_grids.csv", combo, nrow(final_result))
    write.csv(final_result, file.path(output_directory, output_file_name), row.names = FALSE)
    
    cat(sprintf("  -> 已保存: %s\n", output_file_name))
    
  } else {
    missing_files <- bio_files[!file.exists(bio_files)]
    warning(sprintf("情景 %s 缺少文件，跳过处理。\n缺少的第一个文件示例: %s", combo, missing_files[1]))
  }
}
