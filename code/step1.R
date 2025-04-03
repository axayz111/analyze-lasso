# 加载必要的库
library(readr)
library(ggplot2)
library(dplyr)
library(factoextra)
library(tidyr)
library(umap)       # UMAP降维
library(limma)      # 差异表达分析
library(VennDiagram) # Venn图
library(scales)     # 颜色比例
library(spatstat.geom)  # 替代 maptools
library(tidyverse)

# 设置路径
output_dir <- "/Users/zhuzhuxia/Desktop/1"
data_dir <- output_dir
dir.create(data_dir, showWarnings = FALSE)

# 读取数据
train_data <- read_csv("/Users/zhuzhuxia/Desktop/database/11.csv")
train_group <- read_csv("/Users/zhuzhuxia/Desktop/database/12。csv")
train_id_to_gene <- read_csv("/Users/zhuzhuxia/Desktop/database/13.csv")
valid_data <- read_csv("/Users/zhuzhuxia/Desktop/database/21.csv")
valid_group <- read_csv("/Users/zhuzhuxia/Desktop/database/22.csv")
valid_id_to_gene <- read_csv("/Users/zhuzhuxia/Desktop/database/23.csv")


# 1. 基因映射并删除 id 列
map_gene_and_remove_id <- function(data, id_to_gene) {
  data %>%
    left_join(id_to_gene, by = "id") %>% # 根据 id 列映射基因名
    select(gene, everything(), -id) %>% # 将 gene 列放到第一列，并删除 id 列
    filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% # 删除 gene 列为空或包含空白字符的行
    distinct(gene, .keep_all = TRUE) # 去除重复的基因名
}

train_data <- map_gene_and_remove_id(train_data, train_id_to_gene)
valid_data <- map_gene_and_remove_id(valid_data, valid_id_to_gene)

# 2. 检查并转换缺失值
train_data <- train_data %>%
  mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>%
  mutate(across(-gene, as.numeric))

valid_data <- valid_data %>%
  mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>%
  mutate(across(-gene, as.numeric))

# 3. 数据预处理：填充缺失值
preprocess_data <- function(data, file_name) {
  data <- data %>%
    mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) # 填充缺失值

  if (!grepl("log", file_name)) {
    data_numeric <- data %>% select_if(is.numeric)
    data_numeric <- log1p(data_numeric) # log(x + 1) 变换
    data_numeric <- scale(data_numeric) # 标准化
    data <- bind_cols(data %>% select(gene), data_numeric) # 重新组合
  }

  return(data)
}

train_scaled <- preprocess_data(train_data, "train_data.csv")
valid_scaled <- preprocess_data(valid_data, "valid_data.csv")

# 4. 获取交集基因
common_genes <- intersect(train_scaled$gene, valid_scaled$gene)

# 过滤数据，只保留交集基因
train_common <- train_scaled %>% filter(gene %in% common_genes) %>% arrange(gene)
valid_common <- valid_scaled %>% filter(gene %in% common_genes) %>% arrange(gene)

# 检查两个 common 的基因是否完全一致
if (all(train_common$gene == valid_common$gene)) {
  print("两个 common 的基因完全一样")
} else {
  print("两个 common 的基因不完全一样")
  print(paste("train_common 基因数量:", length(train_common$gene)))
  print(paste("valid_common 基因数量:", length(valid_common$gene)))
  print("train_common 中不在 valid_common 的基因:")
  print(setdiff(train_common$gene, valid_common$gene))
  print("valid_common 中不在 train_common 的基因:")
  print(setdiff(valid_common$gene, train_common$gene))
}

# 5. 保存处理后的数据
if (nrow(train_common) > 0) {
  write_csv(train_common, file.path(data_dir, "traindata.csv"))
}

if (nrow(valid_common) > 0) {
  write_csv(valid_common, file.path(data_dir, "validdata.csv"))
}
# ========== 差异表达分析 ==========
visualize_data <- function(data_scaled, group_info, prefix, output_dir) {
  # 构建设计矩阵并拟合线性模型
  design <- model.matrix(~ factor(group_info$Group))
  fit <- eBayes(lmFit(data_scaled %>% select(-gene), design))

  # 获取差异表达基因表格并筛选
  tT2 <- topTable(fit, coef = 2, adjust.method = "BH", number = Inf)
  dT <- decideTests(fit, adjust.method = "BH", p.value = 0.05)

  # 提取表达矩阵和分组信息
  ex <- data_scaled %>% select(-gene)
  gs <- factor(group_info$Group)

  # 统计图表输出
  pdf(file.path(output_dir, paste0(prefix, "_statistical_plots.pdf")), width = 8, height = 6)
  hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj", ylab = "Number of genes", main = "P-adj value distribution")
  vennDiagram(dT, circle.col = palette())
  qqt(fit$t[!is.na(fit$F)], fit$df.total[!is.na(fit$F)], main = "Moderated t statistic")
  volcanoplot(fit, coef = 1, main = colnames(fit)[1], pch = 20, highlight = sum(dT[, 1] != 0), names = rep('+', nrow(fit)))
  plotMD(fit, column = 1, status = dT[, 1], legend = FALSE, pch = 20, cex = 1); abline(h = 0)
  boxplot(ex[, order(gs)], boxwex = 0.6, notch = TRUE, outline = FALSE, las = 2, col = gs[order(gs)])
  legend("topleft", levels(gs), fill = palette(), bty = "n")
  plotDensities(ex, group = gs, main = "Expression Value Distribution", legend = "topright")
  dev.off()

  # UMAP 降维
  pdf(file.path(output_dir, paste0(prefix, "_umap_plot.pdf")), width = 8, height = 6)
  ump <- umap(t(na.omit(ex[!duplicated(ex), ])), n_neighbors = 15, random_state = 123)
  plot(ump$layout, main = "UMAP plot, nbrs=15", xlab = "", ylab = "", col = gs, pch = 20, cex = 1.5)
  legend("topright", inset = c(-0.15, 0), legend = levels(gs), pch = 20, col = 1:nlevels(gs), title = "Group", pt.cex = 1.5)
  dev.off()

  # 均值-方差趋势
  pdf(file.path(output_dir, paste0(prefix, "_mean_variance_trend.pdf")), width = 8, height = 6)
  plotSA(fit, main = "Mean variance trend")
  dev.off()
}

# 对 train 和 valid 进行可视化分析
visualize_data(train_scaled, train_group, "train", output_dir)
visualize_data(valid_scaled, valid_group, "valid", output_dir)

# 加载必要的库
library(ggplot2)
library(ggpubr)
library(umap)

# 1. 处理异常值
remove_outliers <- function(data) {
  data %>%
    mutate(across(where(is.numeric), ~ {
      Q1 <- quantile(., 0.25, na.rm = TRUE)
      Q3 <- quantile(., 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      lower_bound <- Q1 - 1.5 * IQR
      upper_bound <- Q3 + 1.5 * IQR
      ifelse(. < lower_bound | . > upper_bound, NA, .)  # 将异常值替换为 NA
    }))
}

# 2. 生成多个图表的函数
generate_plots <- function(data, group_info, prefix, output_dir) {
  ex <- data %>% select(-gene)  # 移除 gene 列
  gs <- factor(group_info$Group)

  pdf(file.path(output_dir, paste0(prefix, "_all_plots_no_outliers.pdf")), width = 8, height = 6)

  # (1) 箱线图 Boxplot
  boxplot(ex[, order(gs)], boxwex = 0.6, notch = TRUE, outline = FALSE, las = 2, col = gs[order(gs)])
  legend("topleft", levels(gs), fill = palette(), bty = "n")

  # (2) Q-Q 图
  qqnorm(unlist(ex), main = paste(prefix, "Q-Q Plot"))
  qqline(unlist(ex), col = "red")

  # (3) 火山图 Volcano plot
  logFC <- rowMeans(ex[gs == "case", ]) - rowMeans(ex[gs == "control", ])
  pval <- apply(ex, 1, function(x) t.test(x[gs == "case"], x[gs == "control"])$p.value)
  volcano_data <- data.frame(logFC = logFC, pval = -log10(pval))
  ggplot(volcano_data, aes(x = logFC, y = pval)) +
    geom_point(alpha = 0.5) + theme_minimal() +
    ggtitle(paste(prefix, "Volcano Plot"))

  # (4) MD 图（Mean-Difference plot）
  mean_values <- rowMeans(ex)
  md_data <- data.frame(mean = mean_values, diff = logFC)
  ggplot(md_data, aes(x = mean, y = diff)) +
    geom_point(alpha = 0.5) + theme_minimal() +
    ggtitle(paste(prefix, "MD Plot"))

  # (5) 表达值密度分布图
  plot(density(unlist(ex), na.rm = TRUE), main = paste(prefix, "Density Plot"), col = "blue")

  # (6) UMAP 降维图
  umap_result <- umap(t(ex))
  umap_df <- data.frame(UMAP1 = umap_result$layout[,1], UMAP2 = umap_result$layout[,2], Group = gs)
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
    geom_point() + theme_minimal() + ggtitle(paste(prefix, "UMAP Plot"))

  # (7) 均值-方差趋势图
  variance_values <- apply(ex, 1, var, na.rm = TRUE)
  mean_var_data <- data.frame(mean = mean_values, variance = variance_values)
  ggplot(mean_var_data, aes(x = mean, y = variance)) +
    geom_point(alpha = 0.5) + theme_minimal() +
    ggtitle(paste(prefix, "Mean-Variance Trend"))

  dev.off()
}

# 3. 处理 train 和 valid 数据
train_scaled_no_outliers <- remove_outliers(train_scaled)
valid_scaled_no_outliers <- remove_outliers(valid_scaled)

# 4. 检查是否需要重新绘图
if (!all.equal(train_scaled_no_outliers, train_scaled, check.attributes = FALSE)) {
  message("Train 数据发生变化，重新生成所有图表")
  generate_plots(train_scaled_no_outliers, train_group, "train", output_dir)
} else {
  message("Train 数据未发生变化，不生成新的图表")
}

if (!all.equal(valid_scaled_no_outliers, valid_scaled, check.attributes = FALSE)) {
  message("Valid 数据发生变化，重新生成所有图表")
  generate_plots(valid_scaled_no_outliers, valid_group, "valid", output_dir)
} else {
  message("Valid 数据未发生变化，不生成新的图表")
}