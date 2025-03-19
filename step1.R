library(shiny)
library(readr)
library(ggplot2)
library(dplyr)
library(factoextra)
library(tidyr)
library(umap)
library(limma)
library(VennDiagram)
library(scales)
library(spatstat.geom)
library(tidyverse)
library(ggpubr)
library(Cairo)

ui <- fluidPage(
  titlePanel("Data Processing and Visualization"),

  sidebarLayout(
    sidebarPanel(
      h4("Input Files"),
      textInput("train_data_path", "Train Data Path (.csv):", value = "/Users/zhuzhuxia/Desktop/database/心1-61144valid/GSE61144_series_matrix.csv"),
      textInput("train_group_path", "Train Group Path (.csv):", value = "/Users/zhuzhuxia/Desktop/database/心1-61144valid/分组信息.csv"),
      textInput("train_id_to_gene_path", "Train ID to Gene Path (.csv):", value = "/Users/zhuzhuxia/Desktop/database/心1-61144valid/61144ID2Gene.csv"),
      textInput("valid_data_path", "Validation Data Path (.csv):", value = "/Users/zhuzhuxia/Desktop/database/心2-59867/GSE59867_series_matrix.csv"),
      textInput("valid_group_path", "Validation Group Path (.csv):", value = "/Users/zhuzhuxia/Desktop/database/心2-59867/59867分组信息.csv"),
      textInput("valid_id_to_gene_path", "Validation ID to Gene Path (.csv):", value = "/Users/zhuzhuxia/Desktop/database/心2-59867/GSE59867ID2Gene.csv"),

      h4("Output Settings"),
      textInput("output_dir", "Output Directory:", value = "/Users/zhuzhuxia/Desktop/1"),

      h4("PDF Plot Settings"),
      numericInput("pdf_width", "PDF Width (inches):", value = 8, min = 1),
      numericInput("pdf_height", "PDF Height (inches):", value = 6, min = 1),
      textInput("statistical_plots_title", "Statistical Plots Title:", value = "Statistical Plots"),
      textInput("umap_plot_title", "UMAP Plot Title:", value = "UMAP Plot"),
      textInput("mean_variance_title", "Mean Variance Trend Title:", value = "Mean Variance Trend"),
      textInput("all_plots_no_outliers_title", "All Plots (No Outliers) Title:", value = "All Plots (No Outliers)"),
      textInput("boxplot_color", "Boxplot Color:", value = "steelblue"), # 可选的颜色控制
      textInput("volcano_point_color", "Volcano Plot Point Color:", value = "blue"),
      textInput("umap_point_color", "UMAP Plot Point Color:", value = "blue"),
      textInput("mean_var_point_color", "Mean-Variance Point Color:", value = "blue"),

      actionButton("run_processing", "Run Processing")
    ),

    mainPanel(
      h3("Processing Output"),
      verbatimTextOutput("console_output"),
      downloadButton("download_train_common", "Download Train Common Scaled (.csv)"),
      downloadButton("download_valid_common", "Download Valid Common Scaled (.csv)"),
      downloadButton("download_train_statistical_plots", "Download Train Statistical Plots (.pdf)"),
      downloadButton("download_valid_statistical_plots", "Download Valid Statistical Plots (.pdf)"),
      downloadButton("download_train_umap_plot", "Download Train UMAP Plot (.pdf)"),
      downloadButton("download_valid_umap_plot", "Download Valid UMAP Plot (.pdf)"),
      downloadButton("download_train_mean_variance", "Download Train Mean Variance Trend (.pdf)"),
      downloadButton("download_valid_mean_variance", "Download Valid Mean Variance Trend (.pdf)"),
      downloadButton("download_train_all_plots", "Download Train All Plots (No Outliers) (.pdf)"),
      downloadButton("download_valid_all_plots", "Download Valid All Plots (No Outliers) (.pdf)")
    )
  )
)

server <- function(input, output) {
  output$console_output <- renderPrint({
    input$run_processing

    isolate({
      output_dir <- input$output_dir
      data_dir <- output_dir
      dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

      train_data_path <- input$train_data_path
      train_group_path <- input$train_group_path
      train_id_to_gene_path <- input$train_id_to_gene_path
      valid_data_path <- input$valid_data_path
      valid_group_path <- input$valid_group_path
      valid_id_to_gene_path <- input$valid_id_to_gene_path

      pdf_width <- input$pdf_width
      pdf_height <- input$pdf_height
      statistical_plots_title <- input$statistical_plots_title
      umap_plot_title <- input$umap_plot_title
      mean_variance_title <- input$mean_variance_title
      all_plots_no_outliers_title <- input$all_plots_no_outliers_title
      boxplot_color <- input$boxplot_color
      volcano_point_color <- input$volcano_point_color
      umap_point_color <- input$umap_point_color
      mean_var_point_color <- input$mean_var_point_color

      # 读取数据
      train_data <- read_csv(train_data_path)
      train_group <- read_csv(train_group_path)
      train_id_to_gene <- read_csv(train_id_to_gene_path)
      valid_data <- read_csv(valid_data_path)
      valid_group <- read_csv(valid_group_path)
      valid_id_to_gene <- read_csv(valid_id_to_gene_path)

      # 基因映射和ID移除
      map_gene_and_remove_id <- function(data, id_to_gene) {
        data %>% left_join(id_to_gene, by = "id") %>% select(gene, everything(), -id) %>% filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% distinct(gene, .keep_all = TRUE)
      }
      train_data <- map_gene_and_remove_id(train_data, train_id_to_gene)
      valid_data <- map_gene_and_remove_id(valid_data, valid_id_to_gene)

      # 处理 "null" 值并转换为数值
      train_data <- train_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))
      valid_data <- valid_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))

      # 数据预处理函数
      preprocess_data <- function(data, file_name) {
        data <- data %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
        if (!grepl("log", file_name)) {
          data_numeric <- data %>% select_if(is.numeric)
          data_numeric <- log1p(data_numeric)
          data_numeric <- scale(data_numeric)
          data <- bind_cols(data %>% select(gene), data_numeric)
        }
        return(data)
      }
      train_scaled <- preprocess_data(train_data, "train_data.csv")
      valid_scaled <- preprocess_data(valid_data, "valid_data.csv")

      # 查找共同基因
      common_genes <- intersect(train_scaled$gene, valid_scaled$gene)
      train_common <- train_scaled %>% filter(gene %in% common_genes) %>% arrange(gene)
      valid_common <- valid_scaled %>% filter(gene %in% common_genes) %>% arrange(gene)

      # 保存共同基因数据
      if (nrow(train_common) > 0) write_csv(train_common, file.path(data_dir, "train_common_scaled.csv"))
      if (nrow(valid_common) > 0) write_csv(valid_common, file.path(data_dir, "valid_common_scaled.csv"))

      # 可视化函数
      visualize_data_shiny <- function(data_scaled, group_info, prefix, output_dir, pdf_width, pdf_height, title_prefix) {
        design <- model.matrix(~ factor(group_info$Group))
        fit <- eBayes(lmFit(data_scaled %>% select(-gene), design))
        tT2 <- topTable(fit, coef = 2, adjust.method = "BH", number = Inf)
        dT <- decideTests(fit, adjust.method = "BH", p.value = 0.05)
        ex <- data_scaled %>% select(-gene)
        gs <- factor(group_info$Group)
        pdf(file.path(output_dir, paste0(prefix, "_statistical_plots.pdf")), width = pdf_width, height = pdf_height)
        hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj", ylab = "Number of genes", main = paste(title_prefix, "P-adj value distribution"))
        vennDiagram(dT, circle.col = palette())
        qqt(fit$t[!is.na(fit$F)], fit$df.total[!is.na(fit$F)], main = paste(title_prefix, "Moderated t statistic"))
        volcanoplot(fit, coef = 1, main = paste(title_prefix, colnames(fit)[1]), pch = 20, highlight = sum(dT[, 1] != 0), names = rep('+', nrow(fit)))
        plotMD(fit, column = 1, status = dT[, 1], legend = FALSE, pch = 20, cex = 1); abline(h = 0)
        boxplot(ex[, order(gs)], boxwex = 0.6, notch = TRUE, outline = FALSE, las = 2, col = gs[order(gs)])
        legend("topleft", levels(gs), fill = palette(), bty = "n")
        plotDensities(ex, group = gs, main = paste(title_prefix, "Expression Value Distribution"), legend = "topright")
        dev.off()
        pdf(file.path(output_dir, paste0(prefix, "_umap_plot.pdf")), width = pdf_width, height = pdf_height)
        ump <- umap(t(na.omit(ex[!duplicated(ex), ])), n_neighbors = 15, random_state = 123)
        plot(ump$layout, main = paste(title_prefix, "UMAP plot, nbrs=15"), xlab = "", ylab = "", col = gs, pch = 20, cex = 1.5)
        legend("topright", inset = c(-0.15, 0), legend = levels(gs), pch = 20, col = 1:nlevels(gs), title = "Group", pt.cex = 1.5)
        dev.off()
        pdf(file.path(output_dir, paste0(prefix, "_mean_variance_trend.pdf")), width = pdf_width, height = pdf_height)
        plotSA(fit, main = paste(title_prefix, "Mean variance trend"))
        dev.off()
      }

      visualize_data_shiny(train_scaled, train_group, "train", output_dir, pdf_width, pdf_height, statistical_plots_title)
      visualize_data_shiny(valid_scaled, valid_group, "valid", output_dir, pdf_width, pdf_height, statistical_plots_title)

      remove_outliers <- function(data) {
        data %>% mutate(across(where(is.numeric), ~ {
          Q1 <- quantile(., 0.25, na.rm = TRUE)
          Q3 <- quantile(., 0.75, na.rm = TRUE)
          IQR <- Q3 - Q1
          lower_bound <- Q1 - 1.5 * IQR
          upper_bound <- Q3 + 1.5 * IQR
          ifelse(. < lower_bound | . > upper_bound, NA, .)
        }))
      }

      generate_plots_shiny <- function(data, group_info, prefix, output_dir, pdf_width, pdf_height, title_prefix, boxplot_color, volcano_point_color, umap_point_color, mean_var_point_color) {
        ex <- data %>% select(-gene)
        gs <- factor(group_info$Group)
        pdf(file.path(output_dir, paste0(prefix, "_all_plots_no_outliers.pdf")), width = pdf_width, height = pdf_height)
        boxplot(ex[, order(gs)], boxwex = 0.6, notch = TRUE, outline = FALSE, las = 2, col = boxplot_color); legend("topleft", levels(gs), fill = palette(), bty = "n")
        qqnorm(unlist(ex), main = paste(title_prefix, "Q-Q Plot")); qqline(unlist(ex), col = "red")
        logFC <- rowMeans(ex[gs == "case", ]) - rowMeans(ex[gs == "control", ])
        pval <- apply(ex, 1, function(x) t.test(x[gs == "case"], x[gs == "control"])$p.value)
        volcano_data <- data.frame(logFC = logFC, pval = -log10(pval))
        print(ggplot(volcano_data, aes(x = logFC, y = pval)) + geom_point(alpha = 0.5, color = volcano_point_color) + theme_minimal() + ggtitle(paste(title_prefix, "Volcano Plot")))
        mean_values <- rowMeans(ex)
        md_data <- data.frame(mean = mean_values, diff = logFC)
        print(ggplot(md_data, aes(x = mean, y = diff)) + geom_point(alpha = 0.5, color = volcano_point_color) + theme_minimal() + ggtitle(paste(title_prefix, "MD Plot")))
        plot(density(unlist(ex), na.rm = TRUE), main = paste(title_prefix, "Density Plot"), col = "blue")
        umap_result <- umap(t(ex))
        umap_df <- data.frame(UMAP1 = umap_result$layout[,1], UMAP2 = umap_result$layout[,2], Group = gs)
        print(ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) + geom_point(color = umap_point_color) + theme_minimal() + ggtitle(paste(title_prefix, "UMAP Plot")))
        variance_values <- apply(ex, 1, var, na.rm = TRUE)
        mean_var_data <- data.frame(mean = mean_values, variance = variance_values)
        print(ggplot(mean_var_data, aes(x = mean, y = variance)) + geom_point(alpha = 0.5, color = mean_var_point_color) + theme_minimal() + ggtitle(paste(title_prefix, "Mean-Variance Trend")))
        dev.off()
      }

      train_scaled_no_outliers <- remove_outliers(train_scaled)
      valid_scaled_no_outliers <- remove_outliers(valid_scaled)

      generate_plots_shiny(train_scaled_no_outliers, train_group, "train", output_dir, pdf_width, pdf_height, all_plots_no_outliers_title, boxplot_color, volcano_point_color, umap_point_color, mean_var_point_color)
      generate_plots_shiny(valid_scaled_no_outliers, valid_group, "valid", output_dir, pdf_width, pdf_height, all_plots_no_outliers_title, boxplot_color, volcano_point_color, umap_point_color, mean_var_point_color)

      cat("数据处理和可视化完成，结果已保存到目录：", output_dir, "\n")
    })
  })

  output$download_train_common <- downloadHandler(
    filename = function() {
      paste0("train_common_scaled_", Sys.Date(), ".csv")
    },
    content = function(file) {
      train_data_path <- input$train_data_path
      train_group_path <- input$train_group_path
      train_id_to_gene_path <- input$train_id_to_gene_path
      valid_data_path <- input$valid_data_path
      valid_group_path <- input$valid_group_path
      valid_id_to_gene_path <- input$valid_id_to_gene_path

      train_data <- read_csv(train_data_path)
      train_group <- read_csv(train_group_path)
      train_id_to_gene <- read_csv(train_id_to_gene_path)
      valid_data <- read_csv(valid_data_path)
      valid_group <- read_csv(valid_group_path)
      valid_id_to_gene <- read_csv(valid_id_to_gene_path)

      map_gene_and_remove_id <- function(data, id_to_gene) {
        data %>% left_join(id_to_gene, by = "id") %>% select(gene, everything(), -id) %>% filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% distinct(gene, .keep_all = TRUE)
      }
      train_data <- map_gene_and_remove_id(train_data, train_id_to_gene)
      valid_data <- map_gene_and_remove_id(valid_data, valid_id_to_gene)

      train_data <- train_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))
      valid_data <- valid_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))

      preprocess_data <- function(data, file_name) {
        data <- data %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
        if (!grepl("log", file_name)) {
          data_numeric <- data %>% select_if(is.numeric)
          data_numeric <- log1p(data_numeric)
          data_numeric <- scale(data_numeric)
          data <- bind_cols(data %>% select(gene), data_numeric)
        }
        return(data)
      }
      train_scaled <- preprocess_data(train_data, "train_data.csv")
      valid_scaled <- preprocess_data(valid_data, "valid_data.csv")

      common_genes <- intersect(train_scaled$gene, valid_scaled$gene)
      train_common <- train_scaled %>% filter(gene %in% common_genes) %>% arrange(gene)
      write_csv(train_common, file)
    }
  )

  output$download_valid_common <- downloadHandler(
    filename = function() {
      paste0("valid_common_scaled_", Sys.Date(), ".csv")
    },
    content = function(file) {
      train_data_path <- input$train_data_path
      train_group_path <- input$train_group_path
      train_id_to_gene_path <- input$train_id_to_gene_path
      valid_data_path <- input$valid_data_path
      valid_group_path <- input$valid_group_path
      valid_id_to_gene_path <- input$valid_id_to_gene_path

      train_data <- read_csv(train_data_path)
      train_group <- read_csv(train_group_path)
      train_id_to_gene <- read_csv(train_id_to_gene_path)
      valid_data <- read_csv(valid_data_path)
      valid_group <- read_csv(valid_group_path)
      valid_id_to_gene <- read_csv(valid_id_to_gene_path)

      map_gene_and_remove_id <- function(data, id_to_gene) {
        data %>% left_join(id_to_gene, by = "id") %>% select(gene, everything(), -id) %>% filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% distinct(gene, .keep_all = TRUE)
      }
      train_data <- map_gene_and_remove_id(train_data, train_id_to_gene)
      valid_data <- map_gene_and_remove_id(valid_data, valid_id_to_gene)

      train_data <- train_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))
      valid_data <- valid_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))

      preprocess_data <- function(data, file_name) {
        data <- data %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
        if (!grepl("log", file_name)) {
          data_numeric <- data %>% select_if(is.numeric)
          data_numeric <- log1p(data_numeric)
          data_numeric <- scale(data_numeric)
          data <- bind_cols(data %>% select(gene), data_numeric)
        }
        return(data)
      }
      train_scaled <- preprocess_data(train_data, "train_data.csv")
      valid_scaled <- preprocess_data(valid_data, "valid_data.csv")

      common_genes <- intersect(train_scaled$gene, valid_scaled$gene)
      valid_common <- valid_scaled %>% filter(gene %in% common_genes) %>% arrange(gene)
      write_csv(valid_common, file)
    }
  )

  output$download_train_statistical_plots <- downloadHandler(
    filename = function() { paste0("train_statistical_plots_", Sys.Date(), ".pdf") },
    content = function(file) {
      train_data_path <- input$train_data_path
      train_group_path <- input$train_group_path
      train_id_to_gene_path <- input$train_id_to_gene_path

      train_data <- read_csv(train_data_path)
      train_group <- read_csv(train_group_path)
      train_id_to_gene <- read_csv(train_id_to_gene_path)

      map_gene_and_remove_id <- function(data, id_to_gene) {
        data %>% left_join(id_to_gene, by = "id") %>% select(gene, everything(), -id) %>% filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% distinct(gene, .keep_all = TRUE)
      }
      train_data <- map_gene_and_remove_id(train_data, train_id_to_gene)
      train_data <- train_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))
      preprocess_data <- function(data, file_name) {
        data <- data %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
        if (!grepl("log", file_name)) {
          data_numeric <- data %>% select_if(is.numeric)
          data_numeric <- log1p(data_numeric)
          data_numeric <- scale(data_numeric)
          data <- bind_cols(data %>% select(gene), data_numeric)
        }
        return(data)
      }
      train_scaled <- preprocess_data(train_data, "train_data.csv")
      visualize_data_shiny_download(train_scaled, train_group, "train", file, input$pdf_width, input$pdf_height, input$statistical_plots_title)
    }
  )

  output$download_valid_statistical_plots <- downloadHandler(
    filename = function() { paste0("valid_statistical_plots_", Sys.Date(), ".pdf") },
    content = function(file) {
      valid_data_path <- input$valid_data_path
      valid_group_path <- input$valid_group_path
      valid_id_to_gene_path <- input$valid_id_to_gene_path

      valid_data <- read_csv(valid_data_path)
      valid_group <- read_csv(valid_group_path)
      valid_id_to_gene <- read_csv(valid_id_to_gene_path)

      map_gene_and_remove_id <- function(data, id_to_gene) {
        data %>% left_join(id_to_gene, by = "id") %>% select(gene, everything(), -id) %>% filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% distinct(gene, .keep_all = TRUE)
      }
      valid_data <- map_gene_and_remove_id(valid_data, valid_id_to_gene)
      valid_data <- valid_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))
      preprocess_data <- function(data, file_name) {
        data <- data %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
        if (!grepl("log", file_name)) {
          data_numeric <- data %>% select_if(is.numeric)
          data_numeric <- log1p(data_numeric)
          data_numeric <- scale(data_numeric)
          data <- bind_cols(data %>% select(gene), data_numeric)
        }
        return(data)
      }
      valid_scaled <- preprocess_data(valid_data, "valid_data.csv")
      visualize_data_shiny_download(valid_scaled, valid_group, "valid", file, input$pdf_width, input$pdf_height, input$statistical_plots_title)
    }
  )

  output$download_train_umap_plot <- downloadHandler(
    filename = function() { paste0("train_umap_plot_", Sys.Date(), ".pdf") },
    content = function(file) {
      train_data_path <- input$train_data_path
      train_group_path <- input$train_group_path
      train_id_to_gene_path <- input$train_id_to_gene_path

      train_data <- read_csv(train_data_path)
      train_group <- read_csv(train_group_path)
      train_id_to_gene <- read_csv(train_id_to_gene_path)

      map_gene_and_remove_id <- function(data, id_to_gene) {
        data %>% left_join(id_to_gene, by = "id") %>% select(gene, everything(), -id) %>% filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% distinct(gene, .keep_all = TRUE)
      }
      train_data <- map_gene_and_remove_id(train_data, train_id_to_gene)
      train_data <- train_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))
      preprocess_data <- function(data, file_name) {
        data <- data %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
        if (!grepl("log", file_name)) {
          data_numeric <- data %>% select_if(is.numeric)
          data_numeric <- log1p(data_numeric)
          data_numeric <- scale(data_numeric)
          data <- bind_cols(data %>% select(gene), data_numeric)
        }
        return(data)
      }
      train_scaled <- preprocess_data(train_data, "train_data.csv")
      generate_umap_plot_shiny_download(train_scaled, train_group, "train", file, input$pdf_width, input$pdf_height, input$umap_plot_title)
    }
  )

  output$download_valid_umap_plot <- downloadHandler(
    filename = function() { paste0("valid_umap_plot_", Sys.Date(), ".pdf") },
    content = function(file) {
      valid_data_path <- input$valid_data_path
      valid_group_path <- input$valid_group_path
      valid_id_to_gene_path <- input$valid_id_to_gene_path

      valid_data <- read_csv(valid_data_path)
      valid_group <- read_csv(valid_group_path)
      valid_id_to_gene <- read_csv(valid_id_to_gene_path)

      map_gene_and_remove_id <- function(data, id_to_gene) {
        data %>% left_join(id_to_gene, by = "id") %>% select(gene, everything(), -id) %>% filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% distinct(gene, .keep_all = TRUE)
      }
      valid_data <- map_gene_and_remove_id(valid_data, valid_id_to_gene)
      valid_data <- valid_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))
      preprocess_data <- function(data, file_name) {
        data <- data %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
        if (!grepl("log", file_name)) {
          data_numeric <- data %>% select_if(is.numeric)
          data_numeric <- log1p(data_numeric)
          data_numeric <- scale(data_numeric)
          data <- bind_cols(data %>% select(gene), data_numeric)
        }
        return(data)
      }
      valid_scaled <- preprocess_data(valid_data, "valid_data.csv")
      generate_umap_plot_shiny_download(valid_scaled, valid_group, "valid", file, input$pdf_width, input$pdf_height, input$umap_plot_title)
    }
  )

  output$download_train_mean_variance <- downloadHandler(
    filename = function() { paste0("train_mean_variance_", Sys.Date(), ".pdf") },
    content = function(file) {
      train_data_path <- input$train_data_path
      train_group_path <- input$train_group_path
      train_id_to_gene_path <- input$train_id_to_gene_path

      train_data <- read_csv(train_data_path)
      train_group <- read_csv(train_group_path)
      train_id_to_gene <- read_csv(train_id_to_gene_path)

      map_gene_and_remove_id <- function(data, id_to_gene) {
        data %>% left_join(id_to_gene, by = "id") %>% select(gene, everything(), -id) %>% filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% distinct(gene, .keep_all = TRUE)
      }
      train_data <- map_gene_and_remove_id(train_data, train_id_to_gene)
      train_data <- train_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))
      preprocess_data <- function(data, file_name) {
        data <- data %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
        if (!grepl("log", file_name)) {
          data_numeric <- data %>% select_if(is.numeric)
          data_numeric <- log1p(data_numeric)
          data_numeric <- scale(data_numeric)
          data <- bind_cols(data %>% select(gene), data_numeric)
        }
        return(data)
      }
      train_scaled <- preprocess_data(train_data, "train_data.csv")
      generate_mean_variance_shiny_download(train_scaled, train_group, "train", file, input$pdf_width, input$pdf_height, input$mean_variance_title)
    }
  )

  output$download_valid_mean_variance <- downloadHandler(
    filename = function() { paste0("valid_mean_variance_", Sys.Date(), ".pdf") },
    content = function(file) {
      valid_data_path <- input$valid_data_path
      valid_group_path <- input$valid_group_path
      valid_id_to_gene_path <- input$valid_id_to_gene_path

      valid_data <- read_csv(valid_data_path)
      valid_group <- read_csv(valid_group_path)
      valid_id_to_gene <- read_csv(valid_id_to_gene_path)

      map_gene_and_remove_id <- function(data, id_to_gene) {
        data %>% left_join(id_to_gene, by = "id") %>% select(gene, everything(), -id) %>% filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% distinct(gene, .keep_all = TRUE)
      }
      valid_data <- map_gene_and_remove_id(valid_data, valid_id_to_gene)
      valid_data <- valid_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))
      preprocess_data <- function(data, file_name) {
        data <- data %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
        if (!grepl("log", file_name)) {
          data_numeric <- data %>% select_if(is.numeric)
          data_numeric <- log1p(data_numeric)
          data_numeric <- scale(data_numeric)
          data <- bind_cols(data %>% select(gene), data_numeric)
        }
        return(data)
      }
      valid_scaled <- preprocess_data(valid_data, "valid_data.csv")
      generate_mean_variance_shiny_download(valid_scaled, valid_group, "valid", file, input$pdf_width, input$pdf_height, input$mean_variance_title)
    }
  )

  output$download_train_all_plots <- downloadHandler(
    filename = function() { paste0("train_all_plots_no_outliers_", Sys.Date(), ".pdf") },
    content = function(file) {
      train_data_path <- input$train_data_path
      train_group_path <- input$train_group_path
      train_id_to_gene_path <- input$train_id_to_gene_path

      train_data <- read_csv(train_data_path)
      train_group <- read_csv(train_group_path)
      train_id_to_gene <- read_csv(train_id_to_gene_path)

      map_gene_and_remove_id <- function(data, id_to_gene) {
        data %>% left_join(id_to_gene, by = "id") %>% select(gene, everything(), -id) %>% filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% distinct(gene, .keep_all = TRUE)
      }
      train_data <- map_gene_and_remove_id(train_data, train_id_to_gene)
      train_data <- train_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))
      preprocess_data <- function(data, file_name) {
        data <- data %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
        if (!grepl("log", file_name)) {
          data_numeric <- data %>% select_if(is.numeric)
          data_numeric <- log1p(data_numeric)
          data_numeric <- scale(data_numeric)
          data <- bind_cols(data %>% select(gene), data_numeric)
        }
        return(data)
      }
      train_scaled <- preprocess_data(train_data, "train_data.csv")
      train_scaled_no_outliers <- remove_outliers(train_scaled)
      generate_plots_shiny_download(train_scaled_no_outliers, train_group, "train", file, input$pdf_width, input$pdf_height, input$all_plots_no_outliers_title, input$boxplot_color, input$volcano_point_color, input$umap_point_color, input$mean_var_point_color)
    }
  )

  output$download_valid_all_plots <- downloadHandler(
    filename = function() { paste0("valid_all_plots_no_outliers_", Sys.Date(), ".pdf") },
    content = function(file) {
      valid_data_path <- input$valid_data_path
      valid_group_path <- input$valid_group_path
      valid_id_to_gene_path <- input$valid_id_to_gene_path

      valid_data <- read_csv(valid_data_path)
      valid_group <- read_csv(valid_group_path)
      valid_id_to_gene <- read_csv(valid_id_to_gene_path)

      map_gene_and_remove_id <- function(data, id_to_gene) {
        data %>% left_join(id_to_gene, by = "id") %>% select(gene, everything(), -id) %>% filter(!is.na(gene) & gene != "" & !grepl("^\\s*$", gene)) %>% distinct(gene, .keep_all = TRUE)
      }
      valid_data <- map_gene_and_remove_id(valid_data, valid_id_to_gene)
      valid_data <- valid_data %>% mutate(across(everything(), ~ ifelse(. == "null", NA, .))) %>% mutate(across(-gene, as.numeric))
      preprocess_data <- function(data, file_name) {
        data <- data %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
        if (!grepl("log", file_name)) {
          data_numeric <- data %>% select_if(is.numeric)
          data_numeric <- log1p(data_numeric)
          data_numeric <- scale(data_numeric)
          data <- bind_cols(data %>% select(gene), data_numeric)
        }
        return(data)
      }
      valid_scaled <- preprocess_data(valid_data, "valid_data.csv")
      valid_scaled_no_outliers <- remove_outliers(valid_scaled)
      generate_plots_shiny_download(valid_scaled_no_outliers, valid_group, "valid", file, input$pdf_width, input$pdf_height, input$all_plots_no_outliers_title, input$boxplot_color, input$volcano_point_color, input$umap_point_color, input$mean_var_point_color)
    }
  )
}

# Helper function for downloading visualize_data plots
visualize_data_shiny_download <- function(data_scaled, group_info, prefix, file_path, pdf_width, pdf_height, title_prefix) {
  design <- model.matrix(~ factor(group_info$Group))
  fit <- eBayes(lmFit(data_scaled %>% select(-gene), design))
  tT2 <- topTable(fit, coef = 2, adjust.method = "BH", number = Inf)
  dT <- decideTests(fit, adjust.method = "BH", p.value = 0.05)
  ex <- data_scaled %>% select(-gene)
  gs <- factor(group_info$Group)
  pdf(file_path, width = pdf_width, height = pdf_height)
  hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj", ylab = "Number of genes", main = paste(title_prefix, "P-adj value distribution"))
  vennDiagram(dT, circle.col = palette())
  qqt(fit$t[!is.na(fit$F)], fit$df.total[!is.na(fit$F)], main = paste(title_prefix, "Moderated t statistic"))
  volcanoplot(fit, coef = 1, main = paste(title_prefix, colnames(fit)[1]), pch = 20, highlight = sum(dT[, 1] != 0), names = rep('+', nrow(fit)))
  plotMD(fit, column = 1, status = dT[, 1], legend = FALSE, pch = 20, cex = 1); abline(h = 0)
  boxplot(ex[, order(gs)], boxwex = 0.6, notch = TRUE, outline = FALSE, las = 2, col = gs[order(gs)])
  legend("topleft", levels(gs), fill = palette(), bty = "n")
  plotDensities(ex, group = gs, main = paste(title_prefix, "Expression Value Distribution"), legend = "topright")
  dev.off()
}

# Helper function for downloading generate_umap_plot
generate_umap_plot_shiny_download <- function(data, group_info, prefix, file_path, pdf_width, pdf_height, title) {
  ex <- data %>% select(-gene)
  gs <- factor(group_info$Group)
  pdf(file_path, width = pdf_width, height = pdf_height)
  ump <- umap(t(na.omit(ex[!duplicated(ex), ])), n_neighbors = 15, random_state = 123)
  plot(ump$layout, main = title, xlab = "", ylab = "", col = gs, pch = 20, cex = 1.5)
  legend("topright", inset = c(-0.15, 0), legend = levels(gs), pch = 20, col = 1:nlevels(gs), title = "Group", pt.cex = 1.5)
  dev.off()
}

# Helper function for downloading generate_mean_variance
generate_mean_variance_shiny_download <- function(data, group_info, prefix, file_path, pdf_width, pdf_height, title) {
  ex <- data %>% select(-gene)
  mean_values <- rowMeans(ex)
  variance_values <- apply(ex, 1, var, na.rm = TRUE)
  mean_var_data <- data.frame(mean = mean_values, variance = variance_values)
  pdf(file_path, width = pdf_width, height = pdf_height)
  print(ggplot(mean_var_data, aes(x = mean, y = variance)) + geom_point(alpha = 0.5) + theme_minimal() + ggtitle(title))
  dev.off()
}

# Helper function for downloading generate_plots
generate_plots_shiny_download <- function(data, group_info, prefix, file_path, pdf_width, pdf_height, title_prefix, boxplot_color, volcano_point_color, umap_point_color, mean_var_point_color) {
  ex <- data %>% select(-gene)
  gs <- factor(group_info$Group)
  pdf(file_path, width = pdf_width, height = pdf_height)
  boxplot(ex[, order(gs)], boxwex = 0.6, notch = TRUE, outline = FALSE, las = 2, col = boxplot_color); legend("topleft", levels(gs), fill = palette(), bty = "n")
  qqnorm(unlist(ex), main = paste(title_prefix, "Q-Q Plot")); qqline(unlist(ex), col = "red")
  logFC <- rowMeans(ex[gs == "case", ]) - rowMeans(ex[gs == "control", ])
  pval <- apply(ex, 1, function(x) t.test(x[gs == "case"], x[gs == "control"])$p.value)
  volcano_data <- data.frame(logFC = logFC, pval = -log10(pval))
  print(ggplot(volcano_data, aes(x = logFC, y = pval)) + geom_point(alpha = 0.5, color = volcano_point_color) + theme_minimal() + ggtitle(paste(title_prefix, "Volcano Plot")))
  mean_values <- rowMeans(ex)
  md_data <- data.frame(mean = mean_values, diff = logFC)
  print(ggplot(md_data, aes(x = mean, y = diff)) + geom_point(alpha = 0.5, color = volcano_point_color) + theme_minimal() + ggtitle(paste(title_prefix, "MD Plot")))
  plot(density(unlist(ex), na.rm = TRUE), main = paste(title_prefix, "Density Plot"), col = "blue")
  umap_result <- umap(t(ex))
  umap_df <- data.frame(UMAP1 = umap_result$layout[,1], UMAP2 = umap_result$layout[,2], Group = gs)
  print(ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) + geom_point(color = umap_point_color) + theme_minimal() + ggtitle(paste(title_prefix, "UMAP Plot")))
  variance_values <- apply(ex, 1, var, na.rm = TRUE)
  mean_var_data <- data.frame(mean = mean_values, variance = variance_values)
  print(ggplot(mean_var_data, aes(x = mean, y = variance)) + geom_point(alpha = 0.5, color = mean_var_point_color) + theme_minimal() + ggtitle(paste(title_prefix, "Mean-Variance Trend")))
  dev.off()
}

shinyApp(ui = ui, server = server)