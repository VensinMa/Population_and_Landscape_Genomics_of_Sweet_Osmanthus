# ================= 1. 环境准备 =================
packages <- c("ggplot2", "tidyr", "dplyr", "readr", "patchwork")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(patchwork)

# ================= 2. 参数设置 =================

setwd("/home/vensin/Rstudio/RStudio/Population_and_Landscape_Genomics_of_Sweet_Osmanthus/structure/faststructure")
pop_file <- "202samples.pop" 

# 【1. Group 排序】(组内顺序)
target_group_order <- c("LCJ", "DRS", "RX", "JMX", "ZJS", "HYX", "HJLL",
                        "DA", "DST", "ZJP", "EJ", 
                        "ST", "FZA", "GHX", "YZY", "YX",
                        "YK", "LX", "LS", "CP", "CX", "SXK", "SLZ", "SL", 
                        "SFZ", "SK", "XNF", "JD", "QDH", "LQ", "XC", "WYL") 

# 【2. Region 排序】(大区域顺序)
target_region_order <- c("Southwest - Yunnan", "Southwest - Guizhou", "Central", "East")

target_ks <- c(2, 3, 4) 

my_palette <- c("#f58073", "#edcb40", "#b4d66b", "#81b4d6",
                "#8491B4", "#F39B7F", "#91D1C2", "#DC0000", 
                "#7E6148", "#B09C85", "#f58073", "#edcb40")

# ================= 3. 读取并清洗分组信息 =================

if (!file.exists(pop_file)) stop("找不到 pop 文件！")
pop_info <- read.table(pop_file, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
pop_info <- pop_info[, 1:3]
colnames(pop_info) <- c("SampleID", "Group", "Region")
pop_info <- pop_info %>% filter(SampleID != "")

# 应用双重排序
pop_info$Region <- factor(pop_info$Region, levels = target_region_order)
valid_groups <- intersect(target_group_order, unique(pop_info$Group))
remaining_groups <- setdiff(unique(pop_info$Group), valid_groups)
pop_info$Group <- factor(pop_info$Group, levels = c(valid_groups, remaining_groups))

# ================= 4. 定义绘图函数 =================

plot_k <- function(k_val, pop_data, colors) {
  
  meanq_file <- paste0("LD_faststructure_K_", k_val, ".", k_val, ".meanQ")
  if (!file.exists(meanq_file)) return(NULL)
  
  # 读取 Q 矩阵
  q_matrix <- read.table(meanq_file, header = FALSE)
  colnames(q_matrix) <- paste0("K", 1:k_val)
  
  if (nrow(q_matrix) != nrow(pop_data)) stop("样本数不匹配")
  
  # 合并
  df <- cbind(pop_data, q_matrix)
  
  # === 核心排序逻辑 ===
  df_sorted <- df %>%
    arrange(Region, Group, desc(K1)) %>% 
    mutate(plot_id = 1:n()) # 重建绘图 ID
  
  # === 【新增】计算分割线位置 ===
  # 逻辑：找到每个 Region 内，除最后一个 Group 外的所有 Group 的最大 plot_id
  separator_data <- df_sorted %>%
    group_by(Region, Group) %>%
    summarise(max_id = max(plot_id), .groups = "drop") %>%
    group_by(Region) %>%
    # 去掉每个 Region 里最后一个 Group（因为那是分面的边缘，不需要画虚线）
    filter(max_id != max(max_id)) %>% 
    mutate(sep = max_id + 0.5) # 线画在两个柱子中间
  
  # === 计算 X 轴标签位置 ===
  axis_labels <- df_sorted %>%
    group_by(Region, Group) %>%
    summarise(center = mean(plot_id), .groups = "drop") %>%
    arrange(center)
  
  # 转长格式
  df_long <- df_sorted %>%
    pivot_longer(cols = starts_with("K"), names_to = "Ancestry", values_to = "Proportion")
  
  # === 绘图 ===
  p <- ggplot(df_long, aes(x = plot_id, y = Proportion, fill = Ancestry)) +
    geom_bar(stat = "identity", width = 1, linewidth = 0) +
    
    # 【新增】添加白色虚线隔断
    # inherit.aes = FALSE 很重要，防止干扰主图层的映射
    geom_vline(data = separator_data, aes(xintercept = sep), 
               color = "white", linetype = "dashed", linewidth = 0.3, inherit.aes = FALSE) +
    
    # 分面
    facet_grid(~Region, scales = "free_x", space = "free", switch = "x") +
    
    scale_fill_manual(values = colors) +
    
    # X 轴标签
    scale_x_continuous(
      breaks = axis_labels$center,
      labels = axis_labels$Group,
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    
    labs(y = paste0("K=", k_val), x = NULL) +
    
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "transparent", color = NA),
      strip.text = element_text(size = 12, face = "bold"),
      strip.placement = "outside",
      
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, color = "black"),
      axis.ticks.x = element_line(linewidth = 0.3),
      
      axis.title.y = element_text(size = 12, face = "bold", angle = 0, vjust = 0.5, color = "black"),
      panel.grid = element_blank(),
      panel.spacing = unit(0.08, "cm"), # 谱系间隔断宽度
      legend.position = "none"
    )
  
  return(p)
}

# ================= 5. 批量绘图 =================

plot_list <- list()
for (k in target_ks) {
  current_colors <- my_palette[1:k] 
  p <- plot_k(k, pop_info, current_colors)
  if (!is.null(p)) plot_list[[paste0("K", k)]] <- p
}

# 拼图
final_plot <- wrap_plots(plot_list, ncol = 1) & 
  theme(plot.margin = margin(5, 5, 5, 5))
final_plot
# ================= 6. 保存 =================
##out_filename <- paste0("Structure_Plot_Sorted_K", paste(target_ks, collapse = "-"), ".pdf")
##ggsave(out_filename, final_plot, width = 15, height = 3 * length(target_ks))

## print(paste("完成！已保存为:", out_filename))
