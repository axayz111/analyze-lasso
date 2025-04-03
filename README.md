📖 README

🌟 Overview

This project processes gene expression data using a structured workflow, including preprocessing, differential expression analysis, and LASSO regression analysis. The step-by-step approach ensures accurate and reproducible results.

🔍 Workflow

1️⃣ Load Necessary Packages

The script begins by loading all required libraries for data manipulation, visualization, and statistical analysis.

2️⃣ Set Paths & Create Output Directory

Define file paths for input data and output results.

Ensure that necessary directories exist for storing outputs.

3️⃣ Read Data

📂 The script loads:

Training and validation datasets.

Grouping information and gene ID mapping files.

4️⃣ Perform Gene Mapping & Remove ID Columns

Define and call the map_gene_and_remove_id function to map gene IDs and remove unnecessary columns for both datasets.

5️⃣ Handle Missing Values

Replace "null" values in both datasets with NA.

Convert all columns (except the "gene" column) to numeric values.

6️⃣ Data Preprocessing

Define and call preprocess_data to:
✅ Fill missing values
✅ Apply log1p transformation
✅ Normalize the data

7️⃣ Gene Intersection & Data Filtering

🔗 Identify common genes in the training and validation sets and filter datasets accordingly.
✅ Ensure the filtered training and validation sets have identical genes.
💾 Save the processed datasets.

8️⃣ Differential Expression Analysis

🔬 Perform differential expression analysis to identify significantly differentially expressed genes.
📊 Visualize results with plots.

9️⃣ Handle Outliers

Define and apply the remove_outliers function to eliminate extreme values.

🔟 Generate Visualization Plots

Define the generate_plots function to create multiple visualizations for both datasets after outlier removal.

1️⃣1️⃣ LASSO Regression Analysis (First Iteration)

📥 Data Preparation:

Re-load part of the data.

Perform LASSO regression analysis:

Transpose, merge, and group training data.

Prepare training data for the LASSO model.

Align genes, merge groups, and prepare validation data.

🎯 Model Training & Evaluation:

Use cross-validation to determine the best lambda values.

Train the LASSO model and predict validation set responses.

Export predictions, actual values, residuals, and model coefficients.

📈 Visualize actual vs. predicted values and residual plots.

1️⃣2️⃣ LASSO Regression Analysis (Second Iteration with Fixed Lambda)

📥 Re-load the data and perform LASSO regression with fixed lambda values.
📊 The preprocessing steps remain the same, but different file paths are used.
📌 Train the LASSO model, predict validation set responses, and export results.
📈 Generate cross-validation curves, LASSO path plots, and regression coefficient visualizations.

📂 Output Files

📜 This project generates:

Processed training and validation datasets.

Differential expression analysis results and plots.

LASSO regression predictions, residuals, and coefficients.

Various visualization plots.

⚙️ Requirements

🛠️ Before running the script, ensure:

Python/R and necessary libraries are installed.

Data files are formatted correctly and paths are set properly.

💡 Notes

⚠️ Check file paths before running the script.
⚙️ Adjust hyperparameters as needed to optimize model performance.

📬 Contact

For any issues or further clarifications, feel free to reach out! 🎉


📖 README | 使用指南

🌟 Overview | 概述
这个项目主要用于处理基因表达数据，采用了一整套科学的工作流程，包括数据预处理、差异表达分析和 LASSO 回归分析。代码逻辑清晰，确保了结果的准确性和可复现性。

🔍 Workflow | 工作流程
1️⃣ Load Necessary Packages | 加载必要的库
首先，我们会加载所有需要的 Python 或 R 依赖库，以便进行数据处理、可视化和统计分析。

2️⃣ Set Paths & Create Output Directory | 设置路径 & 创建输出目录
- 定义输入数据和输出结果的文件路径。
- 检查并创建必要的目录，以确保结果能正确存储。

3️⃣ Read Data | 读取数据
📂 主要包括：
- 训练集和验证集的数据。
- 组信息和基因 ID 映射文件。

4️⃣ Perform Gene Mapping & Remove ID Columns | 基因映射 & 移除 ID 列
- 定义并调用 `map_gene_and_remove_id` 函数，进行基因 ID 映射，并移除不必要的列。

5️⃣ Handle Missing Values | 处理缺失值
- 将数据集中的 "null" 值替换为 `NA`。
- 除“gene”列外，其余所有列转换为数值类型。

6️⃣ Data Preprocessing | 数据预处理
- 定义并调用 `preprocess_data` 函数，对数据进行：
  ✅ 缺失值填充
  ✅ `log1p` 转换
  ✅ 归一化处理

7️⃣ Gene Intersection & Data Filtering | 交集基因筛选 & 数据过滤
🔗 获取训练集和验证集的共同基因，并基于这些基因筛选数据。
✅ 确保训练集和验证集的基因完全一致。
💾 保存处理后的数据。

8️⃣ Differential Expression Analysis | 差异表达分析
🔬 进行差异表达分析，找出显著的差异基因。
📊 生成数据可视化图表。

9️⃣ Handle Outliers | 处理异常值
- 定义并调用 `remove_outliers` 函数，去除数据中的异常值。

🔟 Generate Visualization Plots | 生成可视化图表
- 定义 `generate_plots` 函数，绘制训练集和验证集的多个图表，以便更直观地分析数据。

1️⃣1️⃣ LASSO Regression Analysis (First Iteration) | LASSO 回归分析 (第一次迭代)
📥 **数据准备：**
- 重新读取部分数据。
- 进行 LASSO 回归分析：
  - 转置、合并、分组训练数据。
  - 准备训练数据。
  - 对齐基因、合并组别，准备验证数据。

🎯 **模型训练 & 评估：**
- 采用交叉验证选择最佳 `lambda` 值。
- 训练 LASSO 模型并预测验证集数据。
- 导出预测值、实际值、残差以及模型系数。
- 📈 生成真实值 vs. 预测值的可视化图表。

1️⃣2️⃣ LASSO Regression Analysis (Second Iteration) | LASSO 回归分析 (第二次迭代)
📥 重新读取数据，并使用固定的 `lambda` 值进行 LASSO 回归。
📊 预处理过程与第一次类似，但使用不同的文件路径。
📌 训练 LASSO 模型，预测验证集数据，导出结果。
📈 绘制交叉验证曲线、LASSO 路径图及回归系数图。

📂 Output Files | 输出文件
📜 主要包括：
- 处理后的训练集和验证集数据。
- 差异表达分析结果及图表。
- LASSO 回归的预测结果、残差、回归系数。
- 其他数据可视化图表。

⚙️ Requirements | 运行环境
🛠️ 运行代码前，请确保：
- 你已安装 **Python/R** 及所需的库。
- **数据文件格式正确**，路径设置无误。

💡 Notes | 注意事项
⚠️ 运行代码前，请先检查文件路径是否正确。
⚙️ 你可以调整超参数，以优化模型性能。

📬 Contact | 联系方式
如果遇到问题或需要帮助，欢迎联系项目成员！🎉

# 📖 README

## 🌟 概要
このプロジェクトは、遺伝子発現データを処理するためのワークフローを提供します。前処理、差次発現解析、および LASSO 回帰分析を含み、正確で再現可能な結果を保証します。

## 🔍 ワークフロー
### 1️⃣ 必要なパッケージの読み込み
データ処理、可視化、および統計分析のためのライブラリをロードします。

### 2️⃣ パスの設定と出力ディレクトリの作成
- 入力データおよび出力結果のファイルパスを定義。
- 結果を正しく保存するためのディレクトリを作成。

### 3️⃣ データの読み込み
📂 以下のデータをロードします：
- 訓練データセットと検証データセット。
- グループ情報と遺伝子 ID マッピングファイル。

### 4️⃣ 遺伝子マッピングと ID 列の削除
- `map_gene_and_remove_id` 関数を定義・呼び出し、遺伝子 ID のマッピングを行い、不要な列を削除。

### 5️⃣ 欠損値の処理
- データセット内の "null" 値を `NA` に置換。
- "gene" 列以外のすべての列を数値型に変換。

### 6️⃣ データの前処理
- `preprocess_data` 関数を定義・呼び出し、以下を実施：
  ✅ 欠損値の補完
  ✅ `log1p` 変換
  ✅ 正規化処理

### 7️⃣ 遺伝子の交差とデータフィルタリング
🔗 訓練セットと検証セットの共通遺伝子を取得し、それに基づいてデータをフィルタリング。
✅ フィルタリング後の訓練セットと検証セットの遺伝子が一致することを確認。
💾 処理済みデータを保存。

### 8️⃣ 差次発現解析
🔬 差次発現解析を実行し、重要な遺伝子を特定。
📊 可視化のためのプロットを生成。

### 9️⃣ 外れ値の処理
- `remove_outliers` 関数を定義し、異常値を除去。

### 🔟 可視化プロットの生成
- `generate_plots` 関数を定義し、外れ値処理後のデータセットを視覚化。

### 1️⃣1️⃣ LASSO 回帰分析（第1回）
📥 **データ準備:**
- データを再読み込み。
- LASSO 回帰分析を実施：
  - 訓練データを転置、結合、グループ化。
  - LASSO モデルのためのデータ準備。
  - 遺伝子の整列、グループの統合、検証データの準備。

🎯 **モデルの学習と評価:**
- 交差検証を用いて最適な `lambda` 値を決定。
- LASSO モデルを学習し、検証セットを予測。
- 予測値、実測値、残差、モデル係数をエクスポート。
- 📈 実測値 vs. 予測値のプロットを生成。

### 1️⃣2️⃣ LASSO 回帰分析（固定 `lambda` を使用した第2回）
📥 データを再読み込みし、固定 `lambda` で LASSO 回帰を実施。
📊 前処理は第1回と同様ですが、異なるファイルパスを使用。
📌 LASSO モデルを学習し、検証セットを予測し、結果をエクスポート。
📈 交差検証曲線、LASSO パスプロット、回帰係数プロットを生成。

## 📂 出力ファイル
📜 本プロジェクトでは以下を生成します：
- 処理済みの訓練データセットと検証データセット。
- 差次発現解析の結果とプロット。
- LASSO 回帰の予測値、残差、回帰係数。
- 各種可視化プロット。

## ⚙️ 必要条件
🛠️ スクリプトを実行する前に、以下を確認してください：
- **Python/R** および必要なライブラリがインストールされていること。
- **データファイルの形式が正しいこと** およびパスが適切に設定されていること。

## 💡 注意事項
⚠️ スクリプトを実行する前に、ファイルパスを確認してください。
⚙️ モデルの性能を最適化するためにハイパーパラメータを調整できます。

## 📬 お問い合わせ
問題が発生した場合や詳細な説明が必要な場合は、お気軽にご連絡ください！ 🎉



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
