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

# 定义指标列表 (区分来源)
# hardy: 来自 --hardy (.hwe文件), site-based
# het:   来自 --het   (.het文件), individual-based
metrics <- c("pi", "tajima", "ho_hardy", "he_hardy", "ho_het", "he_het") 

results_list <- list()

# ================= 3. 定义核心处理函数 =================

# --- 辅助函数1: 解析 .hwe 文件 (Site-based) ---
calc_het_from_hwe <- function(df, type) {
  # type 输入为 "ho_hardy" 或 "he_hardy"
  target_col <- if (type == "ho_hardy") "OBS(HOM1/HET/HOM2)" else "E(HOM1/HET/HOM2)"
  
  if (!target_col %in% colnames(df)) return(NA)
  
  raw_counts <- df[[target_col]]
  
  # 快速拆分
  split_counts <- tstrsplit(raw_counts, "/")
  c1 <- as.numeric(split_counts[[1]])
  c2 <- as.numeric(split_counts[[2]]) # Het
  c3 <- as.numeric(split_counts[[3]])
  
  total <- c1 + c2 + c3
  heterozygosity <- ifelse(total == 0, 0, c2 / total)
  
  return(mean(heterozygosity, na.rm = TRUE))
}

# --- 辅助函数2: 解析 .het 文件 (Individual-based) ---
calc_val_from_het <- function(df, type) {
  # type 输入为 "ho_het" 或 "he_het"
  # .het 文件列: INDV, O(HOM), E(HOM), N_SITES, F
  
  if (!all(c("O(HOM)", "E(HOM)", "N_SITES") %in% colnames(df))) return(NA)
  
  # 移除无效位点数为0的行，防止除以0
  df <- df[df$N_SITES > 0, ]
  if (nrow(df) == 0) return(NA)
  
  if (type == "ho_het") {
    # Ho = (总位点 - 观测纯合) / 总位点
    vals <- (df$N_SITES - df[["O(HOM)"]]) / df$N_SITES
  } else {
    # He = (总位点 - 期望纯合) / 总位点
    vals <- (df$N_SITES - df[["E(HOM)"]]) / df$N_SITES
  }
  
  # 返回群体内所有个体的平均值
  return(mean(vals, na.rm = TRUE))
}

# --- 主处理函数 ---
process_files <- function(dataset, level, metric_type) {
  
  target_dir <- file.path(work_dir, dataset, level)
  
  # 1. 确定文件后缀
  suffix <- ""
  if (metric_type == "pi") {
    suffix <- ".windowed.pi"
  } else if (metric_type == "tajima") {
    suffix <- ".Tajima.D"
  } else if (grepl("_hardy", metric_type)) {
    suffix <- ".hwe"
  } else if (grepl("_het", metric_type)) {
    suffix <- ".het"
  }
  
  # 2. 获取文件
  files <- list.files(target_dir, pattern = paste0(suffix, "$"), full.names = TRUE)
  if (length(files) == 0) return(NULL)
  
  # 3. 循环读取处理
  temp_res <- lapply(files, function(f) {
    # 仅移除后缀，不做名称清洗
    group_name <- sub(suffix, "", basename(f), fixed = TRUE)
    
    df <- tryCatch({
      fread(f, header = TRUE, showProgress = FALSE)
    }, error = function(e) return(NULL))
    
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    mean_val <- NA
    
    # 根据指标类型分发计算逻辑
    if (metric_type == "pi") {
      mean_val <- mean(df[["PI"]], na.rm = TRUE)
    } else if (metric_type == "tajima") {
      mean_val <- mean(df[["TajimaD"]], na.rm = TRUE)
    } else if (grepl("_hardy", metric_type)) {
      mean_val <- calc_het_from_hwe(df, metric_type)
    } else if (grepl("_het", metric_type)) {
      mean_val <- calc_val_from_het(df, metric_type)
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

total_steps <- length(datasets) * length(levels) * length(metrics)
current_step <- 0
counter <- 1

for (d in datasets) {
  for (l in levels) {
    for (m in metrics) {
      current_step <- current_step + 1
      cat(sprintf("[%d/%d] 处理中: %s | %s | %s ...\n", current_step, total_steps, d, l, m))
      
      res <- process_files(d, l, m)
      
      if (!is.null(res)) {
        results_list[[counter]] <- res
        counter <- counter + 1
      }
    }
  }
}

# ================= 5. 数据整合与公式计算 =================
cat("\n正在合并数据并计算 Fis...\n")

if (length(results_list) == 0) stop("未读取到任何有效数据。")

final_long <- do.call(rbind, results_list)

# 映射漂亮的列名 (保持后缀区分)
# 这里的 Metric 已经是带后缀的了 (如 ho_hardy)，我们将其格式化首字母大写
# pi -> pi, tajima -> TajimaD
# ho_hardy -> Ho_hardy, he_het -> He_het
metric_map <- function(m) {
  if (m == "pi") return("pi")
  if (m == "tajima") return("TajimaD")
  # 将 ho_hardy 转换为 Ho_hardy
  parts <- strsplit(m, "_")[[1]]
  paste0(toupper(substr(parts[1], 1, 1)), substr(parts[1], 2, nchar(parts[1])), "_", parts[2])
}
final_long$Pretty_Metric <- sapply(final_long$Metric, metric_map)

# 构造最终列名: Dataset_Metric (例如 ALL_SNP_Ho_hardy)
final_long$Target_Col <- paste(final_long$Dataset, final_long$Pretty_Metric, sep = "_")

# 转宽表
final_wide <- final_long %>%
  select(Level, Group, Target_Col, Mean_Value) %>%
  pivot_wider(names_from = Target_Col, values_from = Mean_Value) %>%
  arrange(Level, Group)

# 计算 Fis = 1 - Ho/He (分别针对 Hardy 和 Het 两种来源)
cat("正在应用公式计算 Fis ...\n")
final_wide <- final_wide %>%
  mutate(
    # Hardy 来源
    ALL_SNP_Fis_hardy = 1 - (ALL_SNP_Ho_hardy / ALL_SNP_He_hardy),
    LD_SNP_Fis_hardy  = 1 - (LD_SNP_Ho_hardy / LD_SNP_He_hardy),
    
    # Het 来源
    ALL_SNP_Fis_het   = 1 - (ALL_SNP_Ho_het / ALL_SNP_He_het),
    LD_SNP_Fis_het    = 1 - (LD_SNP_Ho_het / LD_SNP_He_het)
  )

# 定义列顺序
base_metrics <- c("pi", "TajimaD")
hardy_metrics <- c("Ho_hardy", "He_hardy", "Fis_hardy")
het_metrics <- c("Ho_het", "He_het", "Fis_het")

# 构建期望的列名列表
desired_cols <- c("Level", "Group")
for (d in datasets) {
  # 基础指标
  desired_cols <- c(desired_cols, paste(d, base_metrics, sep="_"))
  # Hardy 指标
  desired_cols <- c(desired_cols, paste(d, hardy_metrics, sep="_"))
  # Het 指标
  desired_cols <- c(desired_cols, paste(d, het_metrics, sep="_"))
}

# 排序并保存
cols_to_keep <- intersect(desired_cols, colnames(final_wide))
final_wide <- final_wide[, cols_to_keep]

output_file <- "Genetic_diversity_summary_DualSource.csv"
write.csv(final_wide, output_file, row.names = FALSE)

cat("====================================================\n")
cat("统计完成！结果已保存至:", output_file, "\n")
print(head(final_wide))
