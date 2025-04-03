library(glmnet)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(Cairo)
library(openxlsx)

train_data_path <- "/Users/zhuzhuxia/Desktop/1--3.csv"
valid_data_path <- "/Users/zhuzhuxia/Desktop/2--4.csv"
train_group_path <- "/Users/zhuzhuxia/Desktop/1和3.csv"
valid_group_path <- "/Users/zhuzhuxia/Desktop/2和4.csv"

output_folder <- "/Users/zhuzhuxia/Desktop/3"

train_data <- read_csv(train_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
valid_data <- read_csv(valid_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
train_group <- read_csv(train_group_path) %>% mutate(Sample = as.character(Sample))
valid_group <- read_csv(valid_group_path) %>% mutate(Sample = as.character(Sample))

train_genes <- unique(as.character(train_data$gene))
train_transposed <- train_data %>% select(-gene) %>% t() %>% as.data.frame() %>% setNames(train_genes) %>% rownames_to_column("Sample")
train_combined <- inner_join(train_transposed, train_group, by = "Sample")
X_train <- as.matrix(train_combined[, train_genes])
y_train <- train_combined$Group

valid_transposed <- valid_data %>% select(-gene) %>% t() %>% as.data.frame() %>% setNames(valid_data$gene) %>% rownames_to_column("Sample")
common_genes <- intersect(train_genes, valid_data$gene)
missing_genes <- setdiff(train_genes, common_genes)
if (length(missing_genes) > 0) {
  valid_aligned <- valid_transposed %>% mutate(across(all_of(missing_genes), ~ 0)) %>% select(Sample, all_of(train_genes))
} else {
  valid_aligned <- valid_transposed %>% select(Sample, all_of(train_genes))
}
valid_combined <- left_join(valid_aligned, valid_group, by = "Sample") %>% filter(!is.na(Group))
X_valid <- as.matrix(valid_combined[, train_genes])
y_valid <- valid_combined$Group

set.seed(123)
cv_model <- cv.glmnet(X_train, y_train, alpha = 1)
best_lambda <- cv_model$lambda.min
cat("交叉验证确定的最佳λ值:", best_lambda, "\n")

fixed_lambda <- 0.1
cat("固定的 λ 值被设置为:", fixed_lambda, "\n")

lasso_model <- glmnet(X_train, y_train, alpha = 1, lambda = fixed_lambda)

predictions <- predict(lasso_model, newx = X_valid)
mse <- mean((predictions - y_valid)^2)
cat("验证集MSE (使用固定 λ 值):", mse, "\n")

results <- data.frame(Sample = valid_combined$Sample, Predicted = predictions[, 1], Actual = y_valid)
write.xlsx(results, file = paste0(output_folder, "固定λ值_预测结果.xlsx"))

cat("训练集特征维度:", dim(X_train), "\n")
cat("训练集目标变量维度:", length(y_train), "\n")
cat("验证集特征维度:", dim(X_valid), "\n")
cat("验证集目标变量维度:", length(y_valid), "\n")

CairoPDF(file = paste0(output_folder, "cv_curve.pdf"), width = 10, height = 6)
plot(cv_model)
dev.off()

CairoPDF(file = paste0(output_folder, "固定λ值_lasso_path_plot.pdf"), width = 10, height = 6)
plot(cv_model$glmnet.fit, xvar = "lambda", label = TRUE)
title("LASSO回归路径图 (固定 λ 值)", line = 2.5)
dev.off()

coef_values <- as.matrix(coef(lasso_model))
coef_values <- as.data.frame(coef_values)
coef_values <- coef_values[coef_values[, 1] != 0, , drop = FALSE]
coef_values <- coef_values[-1, , drop = FALSE]
coef_values$gene <- rownames(coef_values)
colnames(coef_values)[1] <- "Coefficient"

CairoPDF(file = paste0(output_folder, "固定λ值_lasso_coefficients.pdf"), width = 10, height = 6)
print(ggplot(coef_values, aes(x = reorder(gene, -Coefficient), y = Coefficient)) + geom_bar(stat = "identity") + coord_flip() + labs(title = "LASSO回归系数 (固定 λ 值)", x = "基因", y = "回归系数") + theme_minimal())
dev.off()

write.xlsx(coef_values, file = paste0(output_folder, "固定λ值_selected_genes.xlsx"), rowNames = TRUE)