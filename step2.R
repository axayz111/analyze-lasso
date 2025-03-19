library(shiny)
library(glmnet)
library(readxl)
library(openxlsx)
library(dplyr)
library(tibble)
library(ggplot2)
library(readr)
library(Cairo)

ui <- fluidPage(
  titlePanel("LASSO Regression Analysis"),

  sidebarLayout(
    sidebarPanel(
      h4("Input Files"),
      textInput("train_data_path", "Train Data Path (.csv):", value = "/Users/zhuzhuxia/Desktop/1--3.csv"),
      textInput("valid_data_path", "Validation Data Path (.csv):", value = "/Users/zhuzhuxia/Desktop/2--4.csv"),
      textInput("train_group_path", "Train Group Path (.csv):", value = "/Users/zhuzhuxia/Desktop/1和3.csv"),
      textInput("valid_group_path", "Validation Group Path (.csv):", value = "/Users/zhuzhuxia/Desktop/2和4.csv"),

      h4("Output Settings"),
      textInput("output_folder", "Output Folder:", value = "/Users/zhuzhuxia/Desktop/3"),

      h4("PDF Plot Settings"),
      numericInput("pdf_width", "PDF Width (inches):", value = 10, min = 1),
      numericInput("pdf_height", "PDF Height (inches):", value = 6, min = 1),
      textInput("plot_title", "Main Plot Title:", value = "LASSO回归系数 (固定 λ 值)"),
      textInput("path_plot_title", "LASSO Path Plot Title:", value = "LASSO回归路径图 (固定 λ 值)"),
      textInput("coefficients_plot_title", "Coefficients Plot Title:", value = "LASSO回归系数 (固定 λ 值)"),
      textInput("pdf_color", "Bar Color:", value = "steelblue"), # 可选的颜色控制

      h4("LASSO Settings"),
      numericInput("fixed_lambda", "Fixed Lambda Value:", value = 0.1),

      actionButton("run_analysis", "Run Analysis")
    ),

    mainPanel(
      h3("Analysis Output"),
      verbatimTextOutput("console_output"),
      downloadButton("download_results", "Download Prediction Results (.xlsx)"),
      downloadButton("download_coefficients", "Download Model Coefficients (.xlsx)"),
      downloadButton("download_selected_genes", "Download Selected Genes (.xlsx)"),
      downloadButton("download_cv_curve_pdf", "Download CV Curve (.pdf)"),
      downloadButton("download_lasso_path_pdf", "Download LASSO Path Plot (.pdf)"),
      downloadButton("download_coefficients_pdf", "Download Coefficients Plot (.pdf)")
    )
  )
)

server <- function(input, output) {
  output$console_output <- renderPrint({
    input$run_analysis

    # 使用 isolate 来获取按钮点击时的参数值
    isolate({
      train_data_path <- input$train_data_path
      valid_data_path <- input$valid_data_path
      train_group_path <- input$train_group_path
      valid_group_path <- input$valid_group_path
      output_folder <- input$output_folder
      pdf_width <- input$pdf_width
      pdf_height <- input$pdf_height
      plot_title <- input$plot_title
      path_plot_title <- input$path_plot_title
      coefficients_plot_title <- input$coefficients_plot_title
      pdf_color <- input$pdf_color
      fixed_lambda <- input$fixed_lambda

      dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

      # 读取数据
      train_data <- read_csv(train_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      valid_data <- read_csv(valid_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      train_group <- read_csv(train_group_path) %>% mutate(Sample = as.character(Sample))
      valid_group <- read_csv(valid_group_path) %>% mutate(Sample = as.character(Sample))

      # 数据预处理
      train_genes <- unique(as.character(train_data$gene))
      train_transposed <- train_data %>% select(-gene) %>% t() %>% as.data.frame() %>% setNames(train_genes) %>% rownames_to_column("Sample")
      train_combined <- inner_join(train_transposed, train_group, by = "Sample")
      X_train <- as.matrix(train_combined[, train_genes])
      y_train <- train_combined$Group

      valid_transposed <- valid_data %>% select(-gene) %>% t() %>% as.data.frame() %>% setNames(valid_data$gene) %>% rownames_to_column("Sample")
      common_genes <- intersect(train_genes, valid_data$gene)
      missing_genes <- setdiff(train_genes, common_genes)
      if (length(missing_genes) > 0) {
        valid_aligned <- valid_transposed %>%
          mutate(across(all_of(missing_genes), ~ 0)) %>%
          select(Sample, all_of(train_genes))
      } else {
        valid_aligned <- valid_transposed %>%
          select(Sample, all_of(train_genes))
      }
      valid_combined <- left_join(valid_aligned, valid_group, by = "Sample") %>%
        filter(!is.na(Group))
      X_valid <- as.matrix(valid_combined[, train_genes])
      y_valid <- valid_combined$Group

      # 训练模型 - 交叉验证确定最佳lambda
      set.seed(123)
      cv_model <- cv.glmnet(X_train, y_train, alpha = 1)
      best_lambda <- cv_model$lambda.min
      cat("交叉验证确定的最佳λ值:", best_lambda, "\n")

      cat("固定的 λ 值被设置为:", fixed_lambda, "\n")

      # 使用固定的λ值进行LASSO回归拟合
      lasso_model <- glmnet(X_train, y_train, alpha = 1, lambda = fixed_lambda)

      # 预测与评估 (使用固定 λ 值模型)
      predictions <- predict(lasso_model, newx = X_valid)
      mse <- mean((predictions - y_valid)^2)
      cat("验证集MSE (使用固定 λ 值):", mse, "\n")

      # 导出结果 (使用固定 λ 值模型)
      results <- data.frame(Sample = valid_combined$Sample, Predicted = predictions[, 1], Actual = y_valid)
      write.xlsx(results, file = file.path(output_folder, "固定λ值_预测结果.xlsx"))

      # 输出特征和目标变量的维度以调试
      cat("训练集特征维度:", dim(X_train), "\n")
      cat("训练集目标变量维度:", length(y_train), "\n")
      cat("验证集特征维度:", dim(X_valid), "\n")
      cat("验证集目标变量维度:", length(y_valid), "\n")

      # 绘制交叉验证曲线并保存为 PDF
      CairoPDF(file = file.path(output_folder, "cv_curve.pdf"), width = pdf_width, height = pdf_height)
      plot(cv_model)
      dev.off()

      # 绘制LASSO路径图并保存为 PDF
      CairoPDF(file = file.path(output_folder, "固定λ值_lasso_path_plot.pdf"), width = pdf_width, height = pdf_height)
      plot(cv_model$glmnet.fit, xvar = "lambda", label = TRUE)
      title(path_plot_title, line = 2.5)
      dev.off()

      # 提取回归系数并绘制图形 (使用固定 λ 值模型)
      coef_values <- as.matrix(coef(lasso_model))
      coef_values <- as.data.frame(coef_values)
      coef_values <- coef_values[coef_values[, 1] != 0, , drop = FALSE]
      coef_values <- coef_values[-1, , drop = FALSE]

      # 添加基因名，并重命名列
      coef_values$gene <- rownames(coef_values)
      colnames(coef_values)[1] <- "Coefficient"

      # 绘制回归系数图并保存为 PDF
      CairoPDF(file = file.path(output_folder, "固定λ值_lasso_coefficients.pdf"), width = pdf_width, height = pdf_height)
      print(ggplot(coef_values, aes(x = reorder(gene, -Coefficient), y = Coefficient)) +
              geom_bar(stat = "identity", fill = pdf_color) +
              coord_flip() +
              labs(title = coefficients_plot_title, x = "基因", y = "回归系数") +
              theme_minimal())
      dev.off()

      # 导出最终筛选的基因到 XLSX 格式 (使用固定 λ 值模型)
      write.xlsx(coef_values, file = file.path(output_folder, "固定λ值_selected_genes.xlsx"), rowNames = TRUE)

      cat("所有分析完成，结果已保存到目录：", output_folder, "\n")
    })
  })

  # 下载处理程序
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("固定λ值_预测结果_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      train_data_path <- input$train_data_path
      valid_data_path <- input$valid_data_path
      train_group_path <- input$train_group_path
      valid_group_path <- input$valid_group_path

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
        valid_aligned <- valid_transposed %>%
          mutate(across(all_of(missing_genes), ~ 0)) %>%
          select(Sample, all_of(train_genes))
      } else {
        valid_aligned <- valid_transposed %>%
          select(Sample, all_of(train_genes))
      }
      valid_combined <- left_join(valid_aligned, valid_group, by = "Sample") %>%
        filter(!is.na(Group))
      X_valid <- as.matrix(valid_combined[, train_genes])
      y_valid <- valid_combined$Group

      fixed_lambda <- input$fixed_lambda
      lasso_model <- glmnet(X_train, y_train, alpha = 1, lambda = fixed_lambda)
      predictions <- predict(lasso_model, newx = X_valid)
      results <- data.frame(Sample = valid_combined$Sample, Predicted = predictions[, 1], Actual = y_valid)
      write.xlsx(results, file)
    }
  )

  output$download_coefficients <- downloadHandler(
    filename = function() {
      paste0("固定λ值_model_coefficients_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      train_data_path <- input$train_data_path
      valid_data_path <- input$valid_data_path
      train_group_path <- input$train_group_path
      valid_group_path <- input$valid_group_path

      train_data <- read_csv(train_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      valid_data <- read_csv(valid_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      train_group <- read_csv(train_group_path) %>% mutate(Sample = as.character(Sample))
      valid_group <- read_csv(valid_group_path) %>% mutate(Sample = as.character(Sample))

      train_genes <- unique(as.character(train_data$gene))
      train_transposed <- train_data %>% select(-gene) %>% t() %>% as.data.frame() %>% setNames(train_genes) %>% rownames_to_column("Sample")
      train_combined <- inner_join(train_transposed, train_group, by = "Sample")
      X_train <- as.matrix(train_combined[, train_genes])
      y_train <- train_combined$Group

      fixed_lambda <- input$fixed_lambda
      lasso_model <- glmnet(X_train, y_train, alpha = 1, lambda = fixed_lambda)
      coef_df <- as.data.frame(as.matrix(coef(lasso_model)))
      coef_df$Gene <- rownames(coef_df)
      colnames(coef_df) <- c("Coefficient", "Gene")
      write.xlsx(coef_df, file, row.names = TRUE)
    }
  )

  output$download_selected_genes <- downloadHandler(
    filename = function() {
      paste0("固定λ值_selected_genes_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      train_data_path <- input$train_data_path
      valid_data_path <- input$valid_data_path
      train_group_path <- input$train_group_path
      valid_group_path <- input$valid_group_path

      train_data <- read_csv(train_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      valid_data <- read_csv(valid_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      train_group <- read_csv(train_group_path) %>% mutate(Sample = as.character(Sample))
      valid_group <- read_csv(valid_group_path) %>% mutate(Sample = as.character(Sample))

      train_genes <- unique(as.character(train_data$gene))
      train_transposed <- train_data %>% select(-gene) %>% t() %>% as.data.frame() %>% setNames(train_genes) %>% rownames_to_column("Sample")
      train_combined <- inner_join(train_transposed, train_group, by = "Sample")
      X_train <- as.matrix(train_combined[, train_genes])
      y_train <- train_combined$Group

      fixed_lambda <- input$fixed_lambda
      lasso_model <- glmnet(X_train, y_train, alpha = 1, lambda = fixed_lambda)
      coef_values <- as.matrix(coef(lasso_model))
      coef_values <- as.data.frame(coef_values)
      coef_values <- coef_values[coef_values[, 1] != 0, , drop = FALSE]
      coef_values <- coef_values[-1, , drop = FALSE]
      coef_values$gene <- rownames(coef_values)
      colnames(coef_values)[1] <- "Coefficient"
      write.xlsx(coef_values, file, rowNames = TRUE)
    }
  )

  output$download_cv_curve_pdf <- downloadHandler(
    filename = function() {
      paste0("cv_curve_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      train_data_path <- input$train_data_path
      valid_data_path <- input$valid_data_path
      train_group_path <- input$train_group_path
      valid_group_path <- input$valid_group_path

      train_data <- read_csv(train_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      valid_data <- read_csv(valid_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      train_group <- read_csv(train_group_path) %>% mutate(Sample = as.character(Sample))
      valid_group <- read_csv(valid_group_path) %>% mutate(Sample = as.character(Sample))

      train_genes <- unique(as.character(train_data$gene))
      train_transposed <- train_data %>% select(-gene) %>% t() %>% as.data.frame() %>% setNames(train_genes) %>% rownames_to_column("Sample")
      train_combined <- inner_join(train_transposed, train_group, by = "Sample")
      X_train <- as.matrix(train_combined[, train_genes])
      y_train <- train_combined$Group

      set.seed(123)
      cv_model <- cv.glmnet(X_train, y_train, alpha = 1)
      CairoPDF(file = file, width = input$pdf_width, height = input$pdf_height)
      plot(cv_model)
      dev.off()
    }
  )

  output$download_lasso_path_pdf <- downloadHandler(
    filename = function() {
      paste0("lasso_path_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      train_data_path <- input$train_data_path
      valid_data_path <- input$valid_data_path
      train_group_path <- input$train_group_path
      valid_group_path <- input$valid_group_path

      train_data <- read_csv(train_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      valid_data <- read_csv(valid_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      train_group <- read_csv(train_group_path) %>% mutate(Sample = as.character(Sample))
      valid_group <- read_csv(valid_group_path) %>% mutate(Sample = as.character(Sample))

      train_genes <- unique(as.character(train_data$gene))
      train_transposed <- train_data %>% select(-gene) %>% t() %>% as.data.frame() %>% setNames(train_genes) %>% rownames_to_column("Sample")
      train_combined <- inner_join(train_transposed, train_group, by = "Sample")
      X_train <- as.matrix(train_combined[, train_genes])
      y_train <- train_combined$Group

      set.seed(123)
      cv_model <- cv.glmnet(X_train, y_train, alpha = 1)
      CairoPDF(file = file, width = input$pdf_width, height = input$pdf_height)
      plot(cv_model$glmnet.fit, xvar = "lambda", label = TRUE)
      title(input$path_plot_title, line = 2.5)
      dev.off()
    }
  )

  output$download_coefficients_pdf <- downloadHandler(
    filename = function() {
      paste0("lasso_coefficients_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      train_data_path <- input$train_data_path
      valid_data_path <- input$valid_data_path
      train_group_path <- input$train_group_path
      valid_group_path <- input$valid_group_path

      train_data <- read_csv(train_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      valid_data <- read_csv(valid_data_path) %>% rename(gene = 1) %>% distinct(gene, .keep_all = TRUE)
      train_group <- read_csv(train_group_path) %>% mutate(Sample = as.character(Sample))
      valid_group <- read_csv(valid_group_path) %>% mutate(Sample = as.character(Sample))

      train_genes <- unique(as.character(train_data$gene))
      train_transposed <- train_data %>% select(-gene) %>% t() %>% as.data.frame() %>% setNames(train_genes) %>% rownames_to_column("Sample")
      train_combined <- inner_join(train_transposed, train_group, by = "Sample")
      X_train <- as.matrix(train_combined[, train_genes])
      y_train <- train_combined$Group

      fixed_lambda <- input$fixed_lambda
      lasso_model <- glmnet(X_train, y_train, alpha = 1, lambda = fixed_lambda)
      coef_values <- as.matrix(coef(lasso_model))
      coef_values <- as.data.frame(coef_values)
      coef_values <- coef_values[coef_values[, 1] != 0, , drop = FALSE]
      coef_values <- coef_values[-1, , drop = FALSE]
      coef_values$gene <- rownames(coef_values)
      colnames(coef_values)[1] <- "Coefficient"

      CairoPDF(file = file, width = input$pdf_width, height = input$pdf_height)
      print(ggplot(coef_values, aes(x = reorder(gene, -Coefficient), y = Coefficient)) +
              geom_bar(stat = "identity", fill = input$pdf_color) +
              coord_flip() +
              labs(title = input$coefficients_plot_title, x = "基因", y = "回归系数") +
              theme_minimal())
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)