The logic of the code to process the data is as follows


Load the necessary packages.
Set the path and create the output directory.
Read the data, grouping information and gene ID mapping for the training and validation sets.
Define and call the function map_gene_and_remove_id to perform gene mapping and ID column removal for the training and validation sets.
Replaces the “null” values in the training and validation sets with NA and converts all columns except the “gene” column to numeric values.
Define and call the function preprocess_data to perform missing value filling, log1p conversion and normalization on the training and validation sets.
Obtain the intersection of the genes of the training and validation sets and filter the data based on the intersected genes.
Check whether the genes of the filtered training set and validation set are identical.
Save the processed training set and validation set data.
Perform differential expression analysis and visualize the results.
Define function remove_outliers to handle outliers.
Define the function generate_plots to generate multiple plots.
Generate plots for the training and validation sets after handling outliers.
Re-read part of the data and perform LASSO regression analysis.
Transpose, merge and group the training data and prepare the data for training the LASSO model.
Transpose, align genes, merge groups, and prepare data for validation of the LASSO model.
Use cross-validation to determine the best lambda values and train the LASSO model.
Predict and evaluate the validation set using the trained model.
Export the predictions, actuals, residuals, and model coefficients.
Visualize the actual vs. predicted values and residual plots.
In the last part of the code, read the data again and perform LASSO regression analysis using fixed lambda values.
Similar data processing as before is performed for the training and validation sets, but different file paths are used.
Train the LASSO model using fixed lambda values.
Predict and evaluate the validation set and export the results and model coefficients.
Plot cross-validation curves, LASSO path plots, and regression coefficients.




代码处理数据的逻辑如下


加载必要的包。
设置路径并创建输出目录。
读取训练集和验证集的数据、分组信息和基因ID映射。
定义并调用函数 map_gene_and_remove_id 对训练集和验证集进行基因映射和ID列删除。
将训练集和验证集中的 "null" 值替换为 NA，并将除 "gene" 列外的其他列转换为数值型。
定义并调用函数 preprocess_data 对训练集和验证集进行缺失值填充、log1p 转换和标准化。
获取训练集和验证集基因的交集，并根据交集基因过滤数据。
检查过滤后的训练集和验证集的基因是否完全一致。
保存处理后的训练集和验证集数据。
进行差异表达分析，并可视化结果。
定义函数 remove_outliers 处理异常值。
定义函数 generate_plots 生成多个图表。
对处理异常值后的训练集和验证集进行图表生成。
重新读取部分数据，进行 LASSO 回归分析。
对训练数据进行转置、合并分组，并准备训练 LASSO 模型的数据。
对验证数据进行转置、对齐基因、合并分组，并准备验证 LASSO 模型的数据。
使用交叉验证确定最佳 lambda 值，并训练 LASSO 模型。
使用训练好的模型对验证集进行预测和评估。
导出预测结果、实际值、残差以及模型系数。
可视化实际值 vs 预测值和残差图。
最后一部分代码中，再次读取数据，并使用固定的 lambda 值进行 LASSO 回归分析。
对训练集和验证集进行与前面类似的数据处理，但使用了不同的文件路径。
使用固定的 lambda 值训练 LASSO 模型。
对验证集进行预测和评估，并导出结果和模型系数。
绘制交叉验证曲线、LASSO 路径图和回归系数图。
