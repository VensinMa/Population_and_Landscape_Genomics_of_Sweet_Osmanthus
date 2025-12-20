# ================= 1. 环境准备 =================
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

# 定义要处理的数据集、层级
datasets <- c("ALL_SNP", "LD_SNP")
levels <- c("Population", "Lineage", "Species")

# 指标列表 (不包含 fis，因为我们稍后用公式算)
metrics <- c("pi", "tajima", "ho", "he") 

results_list <- list()

# ================= 3. 定义核心处理函数 =================

# 计算 Ho 或 He 的辅助函数 (解析 .hwe 文件)
calc_het_from_hwe <- function(df, type) {
  target_col <- if (type == "ho") "OBS(HOM1/HET/HOM2)" else "E(HOM1/HET/HOM2)"
  if (!target_col %in% colnames(df)) return(NA)
  
  raw_counts <- df[[target_col]]
  
  # 使用 data.table 的 tstrsplit 快速拆分字符串
  split_counts <- tstrsplit(raw_counts, "/")
  
  c1 <- as.numeric(split_counts[[1]])
  c2 <- as.numeric(split_counts[[2]]) # Het
  c3 <- as.numeric(split_counts[[3]])
  
  total <- c1 + c2 + c3
  heterozygosity <- ifelse(total == 0, 0, c2 / total)
  
  return(mean(heterozygosity, na.rm = TRUE))
}

process_files <- function(dataset, level, metric_type) {
  
  target_dir <- file.path(work_dir, dataset, level)
  
  # 定义该指标对应的文件后缀
  suffix <- ""
  if (metric_type == "pi") {
    suffix <- ".windowed.pi"
  } else if (metric_type == "tajima") {
    suffix <- ".Tajima.D"
  } else if (metric_type %in% c("ho", "he")) {
    suffix <- ".hwe"
  }
  
  # 获取文件列表
  # fixed = TRUE 确保点号被视为点号，而不是正则通配符
  files <- list.files(target_dir, pattern = paste0(suffix, "$"), full.names = TRUE)
  
  if (length(files) == 0) return(NULL)
  
  temp_res <- lapply(files, function(f) {
    # === 这里的逻辑已简化：仅移除后缀，不做任何名称修正 ===
    # fixed = TRUE 确保只移除完全匹配的后缀字符串
    group_name <- sub(suffix, "", basename(f), fixed = TRUE)
    
    # 读取数据
    df <- tryCatch({
      fread(f, header = TRUE, showProgress = FALSE)
    }, error = function(e) return(NULL))
    
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    mean_val <- NA
    
    if (metric_type == "pi") {
      mean_val <- mean(df[["PI"]], na.rm = TRUE)
    } else if (metric_type == "tajima") {
      mean_val <- mean(df[["TajimaD"]], na.rm = TRUE)
    } else if (metric_type %in% c("ho", "he")) {
      mean_val <- calc_het_from_hwe(df, metric_type)
    }
    
    return(data.frame(Group = group_name, Mean_Value = mean_val))
  })
  
  combined_df <- do.call(rbind, temp_res)
  if (is.null(combined_df)) return(NULL)
  
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
    for (m in metrics) {
      res <- process_files(d, l, m)
      if (!is.null(res)) {
        results_list[[counter]] <- res
        counter <- counter + 1
      }
    }
  }
}

# ================= 5. 数据整合与公式计算 =================
cat("正在合并并计算 Fis...\n")

if (length(results_list) == 0) {
  stop("未读取到任何有效数据，请检查路径或文件生成情况。")
}

final_long <- do.call(rbind, results_list)

# 映射显示名称
metric_map <- c("pi" = "pi", "tajima" = "TajimaD", "ho" = "Ho", "he" = "He")
final_long$Pretty_Metric <- metric_map[final_long$Metric]

final_long$Target_Col <- paste(final_long$Dataset, final_long$Pretty_Metric, sep = "_")

# 转宽表 (Pivot Wider)
final_wide <- final_long %>%
  select(Level, Group, Target_Col, Mean_Value) %>%
  pivot_wider(names_from = Target_Col, values_from = Mean_Value) %>%
  arrange(Level, Group)

# 计算 Fis (使用 Ho/He 公式)
final_wide <- final_wide %>%
  mutate(
    ALL_SNP_Fis = 1 - (ALL_SNP_Ho / ALL_SNP_He),
    LD_SNP_Fis  = 1 - (LD_SNP_Ho / LD_SNP_He)
  )

# 定义最终列顺序
desired_cols <- c("Level", "Group", 
                  "ALL_SNP_pi", "LD_SNP_pi",
                  "ALL_SNP_TajimaD", "LD_SNP_TajimaD",
                  "ALL_SNP_Ho", "LD_SNP_Ho",
                  "ALL_SNP_He", "LD_SNP_He",
                  "ALL_SNP_Fis", "LD_SNP_Fis")

# 只保留存在的列
cols_to_keep <- intersect(desired_cols, colnames(final_wide))
final_wide <- final_wide[, cols_to_keep]

# ================= 6. 保存结果 =================
output_file <- "Genetic_diversity_summary_Simple.csv"
write.csv(final_wide, output_file, row.names = FALSE)

cat("====================================================\n")
cat("统计完成！\n")
cat("结果已保存至:", file.path(work_dir, output_file), "\n")
print(head(final_wide))
