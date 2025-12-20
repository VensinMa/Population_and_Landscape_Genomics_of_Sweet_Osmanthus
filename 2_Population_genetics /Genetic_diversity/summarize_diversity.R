# ================= 1. 环境准备 =================
# 检查并加载必要的包
packages <- c("tidyverse", "data.table")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(tidyverse)
library(data.table)

# ================= 2. 参数设置 =================
# 设置工作目录
work_dir <- "/home/vensin/workspace/snpcalling_wild/12.population_genetics/Genetic_diversity"
setwd(work_dir)

# 定义要处理的数据集、层级和指标
datasets <- c("ALL_SNP", "LD_SNP")
levels <- c("Population", "Lineage", "Species")
# 这里列出内部识别用的 metric 关键词
metrics <- c("pi", "tajima", "fis", "ho", "he") 

# 定义存储结果的列表
results_list <- list()

# ================= 3. 定义核心处理函数 =================

# 计算 Ho 或 He 的辅助函数 (解析 .hwe 文件)
calc_het_from_hwe <- function(df, type) {
  # VCFtools .hwe 文件的列名通常为: CHR POS OBS(HOM1/HET/HOM2) E(HOM1/HET/HOM2) ChiSq P_HWE
  
  target_col <- if (type == "ho") "OBS(HOM1/HET/HOM2)" else "E(HOM1/HET/HOM2)"
  
  # 检查列是否存在
  if (!target_col %in% colnames(df)) return(NA)
  
  # 提取目标列数据 (格式如 "20/5/10")
  raw_counts <- df[[target_col]]
  
  # 使用 data.table 的 tstrsplit 快速拆分字符串
  # split_counts 是一个 list，包含 3 个向量 (Hom1, Het, Hom2)
  split_counts <- tstrsplit(raw_counts, "/")
  
  # 转换为数值
  c1 <- as.numeric(split_counts[[1]]) # Hom1
  c2 <- as.numeric(split_counts[[2]]) # Het (杂合子数量)
  c3 <- as.numeric(split_counts[[3]]) # Hom2
  
  # 计算总个体数 (对于 Expect 来说是总期望数)
  total <- c1 + c2 + c3
  
  # 计算杂合度 = 杂合子数 / 总数
  # 防止分母为 0
  heterozygosity <- ifelse(total == 0, 0, c2 / total)
  
  # 返回所有位点的平均值
  return(mean(heterozygosity, na.rm = TRUE))
}

process_files <- function(dataset, level, metric_type) {
  
  target_dir <- file.path(work_dir, dataset, level)
  
  # 根据指标类型确定文件名后缀和处理逻辑
  suffix <- ""
  
  if (metric_type == "pi") {
    suffix <- ".windowed.pi"
  } else if (metric_type == "tajima") {
    suffix <- ".Tajima.D"
  } else if (metric_type == "fis") {
    suffix <- ".het"
  } else if (metric_type %in% c("ho", "he")) {
    suffix <- ".hwe"
  }
  
  # 获取文件列表
  files <- list.files(target_dir, pattern = paste0(suffix, "$"), full.names = TRUE)
  
  if (length(files) == 0) {
    # 某些组合可能没有文件 (例如 LD 数据集如果不跑 HWE)，不报错，仅跳过
    return(NULL)
  }
  
  # 循环读取每个文件
  temp_res <- lapply(files, function(f) {
    group_name <- sub(suffix, "", basename(f))
    
    # 读取数据
    df <- tryCatch({
      fread(f, header = TRUE, showProgress = FALSE)
    }, error = function(e) return(NULL))
    
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    mean_val <- NA
    
    # --- 根据指标类型计算均值 ---
    if (metric_type == "pi") {
      mean_val <- mean(df[["PI"]], na.rm = TRUE)
      
    } else if (metric_type == "tajima") {
      mean_val <- mean(df[["TajimaD"]], na.rm = TRUE)
      
    } else if (metric_type == "fis") {
      # .het 文件包含每个个体的 F 值，求平均作为群体 Fis
      mean_val <- mean(df[["F"]], na.rm = TRUE)
      
    } else if (metric_type %in% c("ho", "he")) {
      # 调用辅助函数解析 .hwe 文件
      mean_val <- calc_het_from_hwe(df, metric_type)
    }
    
    return(data.frame(Group = group_name, Mean_Value = mean_val))
  })
  
  # 合并
  combined_df <- do.call(rbind, temp_res)
  if (is.null(combined_df)) return(NULL)
  
  # 添加元数据
  combined_df$Dataset <- dataset
  combined_df$Level <- level
  combined_df$Metric <- metric_type
  
  return(combined_df)
}

# ================= 4. 批量处理 =================
cat("开始统计数据 (这可能需要几分钟，因为要解析大量 HWE 文件)...\n")

counter <- 1
for (d in datasets) {
  for (l in levels) {
    for (m in metrics) {
      cat(sprintf("正在处理: %s - %s - %s\n", d, l, m))
      res <- process_files(d, l, m)
      if (!is.null(res)) {
        results_list[[counter]] <- res
        counter <- counter + 1
      }
    }
  }
}

# ================= 5. 数据整合与格式转换 =================
cat("正在合并并转换表格格式...\n")

if (length(results_list) == 0) {
  stop("未读取到任何有效数据，请检查路径或文件生成情况。")
}

final_long <- do.call(rbind, results_list)

# 构造列名: Dataset + Metric (例如: ALL_SNP_pi, LD_SNP_he)
# 注意大小写调整
metric_map <- c("pi" = "pi", "tajima" = "TajimaD", 
                "fis" = "Fis", "ho" = "Ho", "he" = "He")
final_long$Pretty_Metric <- metric_map[final_long$Metric]

final_long$Target_Col <- paste(final_long$Dataset, final_long$Pretty_Metric, sep = "_")

# Pivot Wider
final_wide <- final_long %>%
  select(Level, Group, Target_Col, Mean_Value) %>%
  pivot_wider(names_from = Target_Col, values_from = Mean_Value) %>%
  arrange(Level, Group)

# 定义理想的列顺序
desired_cols <- c("Level", "Group", 
                  "ALL_SNP_pi", "LD_SNP_pi",
                  "ALL_SNP_TajimaD", "LD_SNP_TajimaD",
                  "ALL_SNP_Ho", "LD_SNP_Ho",
                  "ALL_SNP_He", "LD_SNP_He",
                  "ALL_SNP_Fis", "LD_SNP_Fis")

# 只保留存在的列并排序
cols_to_keep <- intersect(desired_cols, colnames(final_wide))
final_wide <- final_wide[, cols_to_keep]

# ================= 6. 保存结果 =================
output_file <- "Genetic_diversity_summary_AllMetrics.csv"
write.csv(final_wide, output_file, row.names = FALSE)

cat("====================================================\n")
cat("统计完成！\n")
cat("结果已保存至:", file.path(work_dir, output_file), "\n")
print(head(final_wide))
