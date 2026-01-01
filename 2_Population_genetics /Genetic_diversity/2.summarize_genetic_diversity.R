# ==============================================================================
# 脚本名称: summarize_genetic_diversity_GlobalHet_Fixed.R
# 功能: 
#   1. 读取分组计算的 Pi 和 Tajima's D
#   2. 读取全局计算的 .het 文件 (Individual-based)
#   3. 将全局 .het 数据映射回分组，并确保 Ho/He 被正确保留
# ==============================================================================

# ================= 1. 环境与参数 =================
packages <- c("tidyverse", "data.table")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(tidyverse)
library(data.table)

# [!] 请修改此处为你的实际工作目录
work_dir <- "/home/vensin/workspace/snpcalling_wild/12.population_genetics/Genetic_diversity"
setwd(work_dir)

# 分组信息文件路径
pop_file_path <- "/home/vensin/workspace/snpcalling_wild/11.vcftools_filter/snp/202samples.pop"

datasets <- c("ALL_SNP", "LD_SNP")
levels_pi <- c("Population", "Lineage", "Species") 

# ================= 2. 读取元数据 =================
cat(">>> 读取样本分组信息...\n")
if (!file.exists(pop_file_path)) stop("找不到 POP_FILE: ", pop_file_path)

pop_info <- fread(pop_file_path, header = FALSE, col.names = c("INDV", "Population", "Lineage"))
pop_info$Species <- "All_Samples"

# ================= 3. 处理函数定义 =================

# --- 函数 A: 处理全局 .het 文件 ---
process_global_het <- function(dataset) {
  het_file <- file.path(work_dir, dataset, "Global_Stats", "All_Samples.het")
  
  if (!file.exists(het_file)) {
    warning(paste("文件不存在:", het_file))
    return(NULL)
  }
  
  cat(paste("Processing Global Het for:", dataset, "\n"))
  
  df_het <- fread(het_file)
  
  # 计算个体指标
  df_het <- df_het %>%
    filter(N_SITES > 0) %>%
    mutate(
      Ho_ind = (N_SITES - `O(HOM)`) / N_SITES,
      He_ind = (N_SITES - `E(HOM)`) / N_SITES,
      F_ind  = F
    )
  
  merged_df <- df_het %>%
    inner_join(pop_info, by = "INDV")
  
  summary_list <- list()
  target_levels <- c("Population", "Lineage", "Species")
  
  for (lvl in target_levels) {
    summ <- merged_df %>%
      group_by(Group = .data[[lvl]]) %>% 
      summarise(
        # 直接命名为 Ho_het 和 He_het，去掉 Mean_ 前缀，以便后续匹配
        Ho_het      = mean(Ho_ind, na.rm = TRUE),
        He_het      = mean(He_ind, na.rm = TRUE),
        MeanF_het   = mean(F_ind, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        Level = lvl,
        Dataset = dataset,
        Fis_het_Wright = 1 - (Ho_het / He_het)
      )
    
    summary_list[[lvl]] <- summ
  }
  
  return(bind_rows(summary_list))
}

# --- 函数 B: 处理分组的 Pi 和 Tajima's D ---
process_grouped_stats <- function(dataset, level, metric) {
  target_dir <- file.path(work_dir, dataset, level)
  suffix <- if (metric == "pi") ".windowed.pi" else ".Tajima.D"
  files <- list.files(target_dir, pattern = paste0(suffix, "$"), full.names = TRUE)
  
  if (length(files) == 0) return(NULL)
  
  res_list <- lapply(files, function(f) {
    grp <- sub(suffix, "", basename(f), fixed = TRUE)
    d <- fread(f, showProgress = FALSE)
    
    val <- NA
    if (metric == "pi") val <- mean(d$PI, na.rm = TRUE)
    if (metric == "tajima") val <- mean(d$TajimaD, na.rm = TRUE)
    
    data.frame(Group = grp, Mean_Value = val)
  })
  
  df <- do.call(rbind, res_list)
  if (!is.null(df)) {
    df$Dataset <- dataset
    df$Level <- level
    df$Metric <- metric
  }
  return(df)
}

# ================= 4. 执行计算 =================

# 1. 处理 Het
cat(">>> [1/2] Mapping Global Het to Groups...\n")
het_results_list <- list()
for (d in datasets) {
  het_results_list[[d]] <- process_global_het(d)
}
het_final <- bind_rows(het_results_list)

# 2. 处理 Pi/TD
cat(">>> [2/2] Processing Pi and Tajima's D...\n")
div_metrics <- c("pi", "tajima")
div_results_list <- list()
idx <- 1

for (d in datasets) {
  for (l in levels_pi) {
    for (m in div_metrics) {
      res <- process_grouped_stats(d, l, m)
      if (!is.null(res)) {
        div_results_list[[idx]] <- res
        idx <- idx + 1
      }
    }
  }
}
div_long <- do.call(rbind, div_results_list)

# ================= 5. 数据整合 =================
cat(">>> Merging results...\n")

# Pi/TD 转宽表
div_wide <- div_long %>%
  mutate(Target_Col = paste(Dataset, ifelse(Metric=="pi", "pi", "TajimaD"), sep = "_")) %>%
  select(Level, Group, Target_Col, Mean_Value) %>%
  pivot_wider(names_from = Target_Col, values_from = Mean_Value)

# Het 转宽表 (确保列名正确生成)
# 现在的列名是: Ho_het, He_het, MeanF_het, Fis_het_Wright
het_wide <- het_final %>%
  pivot_longer(cols = c(Ho_het, He_het, MeanF_het, Fis_het_Wright), 
               names_to = "Metric", values_to = "Value") %>%
  mutate(Target_Col = paste(Dataset, Metric, sep = "_")) %>%
  select(Level, Group, Target_Col, Value) %>%
  pivot_wider(names_from = Target_Col, values_from = Value)

# 合并
final_df <- full_join(div_wide, het_wide, by = c("Level", "Group")) %>%
  arrange(Level, Group)

# ================= 6. 保存 =================
# 定义输出顺序 (确保名字严格匹配)
desired_cols <- c("Level", "Group")
metrics_order <- c("pi", "TajimaD", "Ho_het", "He_het", "MeanF_het", "Fis_het_Wright")

for (d in datasets) {
  desired_cols <- c(desired_cols, paste(d, metrics_order, sep="_"))
}

# 严格保留存在的列
final_df <- final_df[, intersect(desired_cols, colnames(final_df))]

output_file <- "Genetic_diversity_Summary_GlobalBased.csv"
write.csv(final_df, output_file, row.names = FALSE)

cat("====================================================\n")
cat("统计完成！结果已保存至:", output_file, "\n")
print(head(final_df))
