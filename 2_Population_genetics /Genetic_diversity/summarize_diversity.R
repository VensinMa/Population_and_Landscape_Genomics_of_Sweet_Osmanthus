# ================= 1. 环境准备 =================
# 检查并加载必要的包
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)
library(data.table)

# ================= 2. 参数设置 =================
# 设置工作目录 (即包含 ALL_SNP 和 LD_SNP 文件夹的那个目录)
work_dir <- "/home/vensin/workspace/snpcalling_wild/12.population_genetics/Genetic_diversity"
setwd(work_dir)

# 定义要处理的数据集和层级
datasets <- c("ALL_SNP", "LD_SNP")
levels <- c("Population", "Lineage", "Species")

# 定义存储结果的列表
results_list <- list()

# ================= 3. 定义读取函数 =================
# 这个函数会读取指定文件夹下所有的 pi 或 Tajima.D 文件并计算均值
process_files <- function(dataset, level, metric_type) {
  
  # 构建路径: e.g., ALL_SNP/Population
  target_dir <- file.path(work_dir, dataset, level)
  
  # 确定文件后缀和列名
  if (metric_type == "pi") {
    suffix <- ".windowed.pi"
    col_name <- "PI"
  } else {
    suffix <- ".Tajima.D"
    col_name <- "TajimaD"
  }
  
  # 获取文件列表
  files <- list.files(target_dir, pattern = paste0(suffix, "$"), full.names = TRUE)
  
  if (length(files) == 0) {
    warning(paste("在", target_dir, "未找到", metric_type, "文件"))
    return(NULL)
  }
  
  # 循环读取每个文件
  temp_res <- lapply(files, function(f) {
    # 提取组名 (去除路径和后缀)
    # basename: /path/to/CP.windowed.pi -> CP.windowed.pi
    # sub: CP.windowed.pi -> CP
    group_name <- sub(suffix, "", basename(f))
    
    # 读取数据 (使用 fread 加速)
    # quiet = TRUE 防止打印过多信息
    df <- tryCatch({
      fread(f, header = TRUE, showProgress = FALSE)
    }, error = function(e) return(NULL))
    
    if (is.null(df) || nrow(df) == 0) {
      return(data.frame(Group = group_name, Mean_Value = NA))
    }
    
    # 计算均值 (去除 NA 值，TajimaD 经常会有 NaN)
    # 注意：如果文件里全是 NaN，mean 会得到 NaN
    mean_val <- mean(df[[col_name]], na.rm = TRUE)
    
    return(data.frame(Group = group_name, Mean_Value = mean_val))
  })
  
  # 合并列表为一个数据框
  combined_df <- do.call(rbind, temp_res)
  
  # 添加元数据列，方便后续 Pivot
  combined_df$Dataset <- dataset
  combined_df$Level <- level
  combined_df$Metric <- metric_type
  
  return(combined_df)
}

# ================= 4. 批量处理 =================
cat("开始统计数据...\n")

counter <- 1
for (d in datasets) {
  for (l in levels) {
    # 1. 处理 Pi
    cat(sprintf("正在处理: %s - %s - Pi\n", d, l))
    res_pi <- process_files(d, l, "pi")
    if (!is.null(res_pi)) results_list[[counter]] <- res_pi
    counter <- counter + 1
    
    # 2. 处理 Tajima's D
    cat(sprintf("正在处理: %s - %s - Tajima's D\n", d, l))
    res_td <- process_files(d, l, "tajima")
    if (!is.null(res_td)) results_list[[counter]] <- res_td
    counter <- counter + 1
  }
}

# ================= 5. 数据整合与格式转换 =================
cat("正在合并并转换表格格式...\n")

# 合并所有结果
final_long <- do.call(rbind, results_list)

# 此时的数据是长格式 (Long Format)，我们需要把它变宽 (Pivot Wider)
# 目标列名: ALL_SNP_pi, LD_SNP_pi, ALL_SNP_TajimaD, LD_SNP_TajimaD

# 1. 构造新的列名 (Dataset + Metric)
final_long$Target_Col <- paste(final_long$Dataset, 
                               ifelse(final_long$Metric == "pi", "pi", "TajimaD"), 
                               sep = "_")

# 2. 执行透视 (Pivot)
final_wide <- final_long %>%
  select(Level, Group, Target_Col, Mean_Value) %>%
  pivot_wider(names_from = Target_Col, values_from = Mean_Value) %>%
  arrange(Level, Group) # 按层级和组名排序

# 调整列顺序 (Level, Group, ALL_SNP_pi, LD_SNP_pi, ALL_SNP_TajimaD, LD_SNP_TajimaD)
desired_order <- c("Level", "Group", 
                   "ALL_SNP_pi", "LD_SNP_pi", 
                   "ALL_SNP_TajimaD", "LD_SNP_TajimaD")

# 只保留存在的列 (防止某些数据没跑完导致列缺失报错)
existing_cols <- intersect(desired_order, colnames(final_wide))
final_wide <- final_wide[, existing_cols]

# ================= 6. 保存结果 =================
output_file <- "Genetic_diversity_summary.csv"
write.csv(final_wide, output_file, row.names = FALSE)

cat("====================================================\n")
cat("统计完成！\n")
cat("结果已保存至:", file.path(work_dir, output_file), "\n")
print(head(final_wide))
