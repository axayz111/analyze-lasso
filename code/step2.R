# 安装和加载必要的包
# install.packages(c("glmnet", "readxl", "openxlsx", "dplyr", "tibble", "ggplot2"))
library(glmnet)    # 用于LASSO回归
library(readxl)    # 用于读取Excel文件
library(openxlsx)  # 用于写入Excel文件
library(dplyr)     # 用于数据处理
library(tibble)    # 用于数据框操作
library(ggplot2)   # 用于可视化
library(readr)

# 设置文件路径
data_dir <- "/Users/zhuzhuxia/Desktop/2"  # 结果保存目录
dir.create(data_dir, showWarnings = FALSE)  # 创建结果目录（如果不存在）

# 定义文件路径
train_data_path <- "/Users/zhuzhuxia/Desktop/1--3.csv"
valid_data_path <- "/Users/zhuzhuxia/Desktop/2--4.csv"
train_group_path <- "/Users/zhuzhuxia/Desktop/1和3.csv"
valid_group_path <- "/Users/zhuzhuxia/Desktop/2和4.csv"

# 读取数据并剔除重复的基因名
train_data <- read_csv(train_data_path)

# 获取train_data的列数
num_cols_train <- ncol(train_data)

# 将train_data的最后一列重命名为 'gene'
train_data <- train_data %>%
  rename(gene = !!names(.)[num_cols_train]) %>% # 使用 !!names(.)[num_cols_train] 获取最后一列的名称并重命名
  distinct(gene, .keep_all = TRUE)  # 保留唯一的基因名

valid_data <- read_csv(valid_data_path)

# 获取valid_data的列数
num_cols_valid <- ncol(valid_data)

# 将valid_data的最后一列重命名为 'gene'
valid_data <- valid_data %>%
  rename(gene = !!names(.)[num_cols_valid]) %>% # 使用 !!names(.)[num_cols_valid] 获取最后一列的名称并重命名
  distinct(gene, .keep_all = TRUE)  # 保留唯一的基因名


train_group <- read_csv(train_group_path) %>%
  mutate(Sample = as.character(Sample))

valid_group <- read_csv(valid_group_path) %>%
  mutate(Sample = as.character(Sample))

# 转置训练数据
train_genes <- train_data$gene
train_transposed <- train_data %>%
  select(-gene) %>%
  t() %>%
  as.data.frame() %>%
  setNames(train_genes) %>%
  rownames_to_column("Sample")

# 合并训练分组
train_combined <- inner_join(train_transposed, train_group, by = "Sample")
X_train <- as.matrix(train_combined[, train_genes])
y_train <- train_combined$Group

# 处理验证数据
valid_transposed <- valid_data %>%
  select(-gene) %>%
  t() %>%
  as.data.frame() %>%
  setNames(valid_data$gene) %>%
  rownames_to_column("Sample")

common_genes <- intersect(train_genes, valid_data$gene)
valid_aligned <- valid_transposed %>%
  mutate(across(all_of(setdiff(train_genes, common_genes)), ~ 0)) %>%
  select(Sample, all_of(train_genes))

valid_combined <- left_join(valid_aligned, valid_group, by = "Sample") %>%
  filter(!is.na(Group))
X_valid <- as.matrix(valid_combined[, train_genes])
y_valid <- valid_combined$Group

# 训练模型
set.seed(123)
cv_model <- cv.glmnet(X_train, y_train, alpha = 1)
best_lambda <- cv_model$lambda.min
lasso_model <- glmnet(X_train, y_train, alpha = 1, lambda = best_lambda)

# 预测与评估
predictions <- predict(lasso_model, s = best_lambda, newx = X_valid)
mse <- mean((predictions - y_valid) ^ 2)
cat("验证集的均方误差 (MSE):", mse, "\n")

# 导出详细结果
results <- data.frame(Sample = valid_combined$Sample,
                      Predicted = predictions[,1],
                      Actual = y_valid,
                      Residuals = predictions[,1] - y_valid)
write.xlsx(results, file.path(data_dir, "prediction_results.xlsx"))

# 可视化：实际值 vs 预测值
p1 <- ggplot(results, aes(x = Actual, y = Predicted)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "实际值 vs 预测值", x = "实际值", y = "预测值") +
  theme_minimal()
ggsave(file.path(data_dir, "actual_vs_predicted.pdf"), plot = p1)

# 可视化：残差图
p2 <- ggplot(results, aes(x = Sample, y = Residuals)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "残差图", x = "样本", y = "残差") +
  theme_minimal()
ggsave(file.path(data_dir, "residuals_plot.pdf"), plot = p2)

# 输出系数
coef_matrix <- as.matrix(coef(lasso_model))  # 将系数转换为常规矩阵
coef_df <- as.data.frame(coef_matrix)  # 转换为数据框
coef_df$Gene <- rownames(coef_df)  # 添加基因名列
colnames(coef_df) <- c("Coefficient", "Gene")  # 重命名列

# 导出模型系数
write.xlsx(coef_df, file.path(data_dir, "model_coefficients.xlsx"), row.names = TRUE)

cat("所有分析完成，结果已保存到目录：", data_dir, "\n")