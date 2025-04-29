# Step 1: Install Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# ✅ Install devtools only if not already installed
if (!require(devtools, quietly = TRUE)) install.packages("devtools")

# ✅ Install the NOLAS package from GitLab
devtools::install_git("https://gitlab.com/algoromics/nolas.git")
BiocManager::install("curatedTCGAData")

install.packages(c("mixOmics", "caret", "MLmetrics", "pROC", "PRROC", "enrichR"))
# Step 2: Load Required Libraries


library(curatedTCGAData)
library(SummarizedExperiment)
library(nolas)
library(mixOmics)
library(caret)
library(MLmetrics)
library(pROC)
library(PRROC)
library(enrichR)
library(DESeq2)
# Step 3: Download TCGA-BRCA Omics Data (including miRNA)
brca_data <- curatedTCGAData(
  diseaseCode = "BRCA",
  assays = c("RNASeq2GeneNorm", "RPPAArray", "miRNASeqGene"),
  version = "1.1.38",
  dry.run = FALSE
)

# 🧬 Step: Extract assay matrices from BRCA dataset
print(brca_data)

rna_data   <- assays(brca_data[["BRCA_RNASeq2GeneNorm-20160128"]])[[1]]
rppa_data  <- assays(brca_data[["BRCA_RPPAArray-20160128"]])[[1]]
mirna_data <- assays(brca_data[["BRCA_miRNASeqGene-20160128"]])[[1]]

RNASeq_matrix <- as.matrix(rna_data)
RPPA_matrix   <- as.matrix(rppa_data)
miRNA_matrix  <- as.matrix(mirna_data)

# 🧬 Step: Define TCGA barcode extraction function
extract_barcode <- function(x) sapply(strsplit(x, "-"), function(parts) paste(parts[1:3], collapse = "-"))

# 🧬 Step: Extract standardized sample barcodes
samples_rna_base   <- extract_barcode(colnames(RNASeq_matrix))
samples_rppa_base  <- extract_barcode(colnames(RPPA_matrix))
samples_mirna_base <- extract_barcode(colnames(miRNA_matrix))

# 🧬 Step: Find common samples across all three omics layers
common_samples <- Reduce(intersect, list(samples_rna_base, samples_rppa_base, samples_mirna_base))
common_samples <- sort(common_samples)




# Extract clinical metadata
clinical_meta <- colData(brca_data)
samples_clinical_base <- extract_barcode(rownames(clinical_meta))

# Extract raw vital status
labels_raw <- clinical_meta$vital_status
# 🧬 Step: Extract survival labels for common samples
labels <- labels_raw[match(common_samples, samples_clinical_base)]
labels[labels == "0"] <- "Alive"
labels[labels == "1"] <- "Deceased"
labels <- factor(labels)
names(labels) <- common_samples

# 🧬 Step: Subset matrices to common samples only
RNASeq_matrix_common <- RNASeq_matrix[, samples_rna_base %in% common_samples, drop = FALSE]
RPPA_matrix_common   <- RPPA_matrix[, samples_rppa_base %in% common_samples, drop = FALSE]
miRNA_matrix_common  <- miRNA_matrix[, samples_mirna_base %in% common_samples, drop = FALSE]

# 🧬 Step: Standardize column names using TCGA barcodes
colnames(RNASeq_matrix_common) <- extract_barcode(colnames(RNASeq_matrix_common))
colnames(RPPA_matrix_common)   <- extract_barcode(colnames(RPPA_matrix_common))
colnames(miRNA_matrix_common)  <- extract_barcode(colnames(miRNA_matrix_common))
# Combine all three omics blocks
multi_omics_data <- list(
  RNASeq = RNASeq_matrix_common,
  RPPA   = RPPA_matrix_common,
  miRNA  = miRNA_matrix_common
)


#Step 10: Stratified Split (50/50) into Training and Testing Sets
set.seed(123)
train_index   <- createDataPartition(labels, p = 0.5, list = FALSE)
train_samples <- common_samples[train_index]
test_samples  <- common_samples[-train_index]

multi_omics_train <- lapply(multi_omics_data, function(x) x[, train_samples])
multi_omics_test  <- lapply(multi_omics_data, function(x) x[, test_samples])

labels_train <- labels[train_samples]
labels_test  <- labels[test_samples]

# 🧬 Step 11: Apply Variance Stabilizing Transformation (VST) to RNA-Seq Training Data
dds <- DESeqDataSetFromMatrix(
  countData = round(multi_omics_train$RNASeq),
  colData = data.frame(row.names = colnames(multi_omics_train$RNASeq)),
  design = ~1
)
vsd <- vst(dds, blind = TRUE)
rna_train_vst <- t(assay(vsd))
# Convert miRNA-Seq training set to DESeq2 object
dds_mirna <- DESeqDataSetFromMatrix(
  countData = round(multi_omics_train$miRNA),
  colData = data.frame(row.names = colnames(multi_omics_train$miRNA)),
  design = ~1
)
vsd_mirna <- varianceStabilizingTransformation(dds_mirna, blind = TRUE)
mirna_train_vst <- t(assay(vsd_mirna))  # Transpose: samples in rows
mirna_train_centered <- scale(mirna_train_vst, center = TRUE, scale = FALSE)
#Step 12: Mean-Center the Transformed RNA-Seq Data and mirna
# RNA-Seq 
rna_train_centered <- scale(rna_train_vst, center = TRUE, scale = FALSE)
# miRNA-Seq 
mirna_train_centered <- scale(mirna_train_vst, center = TRUE, scale = FALSE)


# 🧬 Step 13: Preprocess the RNA-Seq Test Set and mirna
rna_test_raw <- t(multi_omics_test$RNASeq)[, colnames(rna_train_centered)]
rna_test <- scale(rna_test_raw, center = TRUE, scale = FALSE)
rownames(rna_test) <- colnames(multi_omics_test$RNASeq)
rownames(rna_train_centered) <- colnames(vsd)

mirna_test_raw <- t(multi_omics_test$miRNA)[, colnames(mirna_train_centered)]
mirna_test <- scale(mirna_test_raw, center = TRUE, scale = FALSE)
rownames(mirna_test) <- colnames(multi_omics_test$miRNA)
rownames(mirna_train_centered) <- colnames(multi_omics_train$miRNA)


# 🧬 Step 14: Preprocess RPPA Training Data (Impute NA & Inf)
rppa_train <- t(multi_omics_train$RPPA)
X_rppa <- apply(rppa_train, 2, function(col) {
  col[is.na(col) | is.infinite(col)] <- mean(col, na.rm = TRUE)
  return(col)
})
rownames(X_rppa) <- rownames(rppa_train)

# Step 15: Preprocess RPPA Test Set (Same Cleaning)

rppa_test <- t(multi_omics_test$RPPA)
rppa_test <- apply(rppa_test, 2, function(col) {
  col[is.na(col) | is.infinite(col)] <- mean(col, na.rm = TRUE)
  return(col)
})
rownames(rppa_test) <- colnames(multi_omics_test$RPPA)
rppa_train_matrix <- X_rppa
# 🧬 Step 16: Build Final Multi-Omics Training Block
# 💬 Description:
# Select only those samples that are common across RNA-Seq, RPPA, and the label vector.
# This ensures consistent and matched training input for multi-omics integration.

common_train_samples <- Reduce(intersect, list(
  rownames(rna_train_centered),
  rownames(X_rppa),
  rownames(mirna_train_centered),
  names(labels_train)
))

omics_train <- list(
  RNASeq = rna_train_centered[common_train_samples, ],
  RPPA   = rppa_train_matrix[common_train_samples, ],
  miRNA = mirna_train_centered[common_train_samples, ]
)
survival_train_labels <- labels_train
training_metadata <- data.frame(Survival = factor(survival_train_labels[common_train_samples]))

# 🧬 Step 17: Build Final Multi-Omics Test Block
X_test <- list(
  RNASeq = rna_test,
  RPPA   = rppa_test,
  miRNA  = mirna_test
)
# 🧬 Step 18: Train Final NOLAS Model
# 🧠 Train final NOLAS model with 3 omics layers
res_nolas_final <- nolas(
  X = omics_train,                  # Includes RNASeq, RPPA, and miRNA
  Y = training_metadata,           # Survival labels
  nlv = 2,                         # Number of latent variables
  imputeNA = TRUE,                # Impute missing values
  center = TRUE,                  # Mean-center each block
  scale = FALSE,                  # No scaling
  response = "Survival",         # Classification task
  rvc = 0.25                      # Relative variance cutoff for feature selection
)
cat("✅ Final NOLAS model trained on", nrow(omics_train$RNASeq), "samples with", length(omics_train), "omics layers\n")

# 🧪 Step 19: Align Test Blocks with Training Features
# ✅ Define helper functions for test block cleaning and alignment
clean_block <- function(block) {
  block <- as.matrix(block)
  block[!is.finite(block)] <- NA
  
  # Remove columns with all NA values
  block <- block[, colSums(is.na(block)) < nrow(block), drop = FALSE]
  
  # Impute missing values with column mean
  for (j in seq_len(ncol(block))) {
    if (anyNA(block[, j])) {
      mean_value <- mean(block[, j], na.rm = TRUE)
      block[, j][is.na(block[, j])] <- ifelse(is.nan(mean_value), 0, mean_value)
    }
  }
  return(block)
}

fix_block_strict <- function(test_block, train_block) {
  common_cols <- intersect(colnames(test_block), colnames(train_block))
  test_block_clean <- test_block[, common_cols, drop = FALSE]
  
  # Reorder columns to match training block
  test_block_clean <- test_block_clean[, match(common_cols, colnames(train_block)), drop = FALSE]
  return(test_block_clean)
}

# 🧪 Step: Clean and align all test blocks including miRNA
omics_test_clean <- list(
  RNASeq = fix_block_strict(clean_block(rna_test), omics_train$RNASeq),
  RPPA   = fix_block_strict(clean_block(rppa_test), omics_train$RPPA),
  miRNA  = fix_block_strict(clean_block(mirna_test), omics_train$miRNA)
)


# 🔮 Step 21: Make predictions using the trained NOLAS model

predictions <- predict(res_nolas_final, omics_test_clean)
str(predictions)
# 🧠 Step 22: Extract predicted labels from the output

y_prob <- predictions$probY
y_pred <- apply(y_prob, 1, function(probs) {
  if (any(is.na(probs))) return(NA)
  colnames(y_prob)[which.max(probs)]
})
y_pred <- factor(y_pred, levels = c("Alive", "Deceased"))
# 🎯 Step 23: Match predicted and true labels
survival_test_labels <- labels_test
y_true <- labels_test[names(y_pred)]

# ✅ Step 24: Filter valid predictions
valid_idx <- !is.na(y_pred) & !is.na(y_true)
y_pred_clean <- y_pred[valid_idx]
y_true_clean <- y_true[valid_idx]
y_score_nolas <- y_prob[valid_idx, "Deceased"]

# 📊 Step 24: View distribution of classes and example scores
cat("📊 Class distribution in filtered test set:\n")
print(table(y_true_clean))

cat("🔍 Example of predicted 'Deceased' scores:\n")
print(head(y_score_nolas))

# 📈 Step 25: Compute AUC for NOLAS predictions
library(pROC)
roc_obj_nolas <- roc(y_true_clean, y_score_nolas, levels = c("Alive", "Deceased"), direction = "<")
auc_val_nolas <- auc(roc_obj_nolas)

library(MLmetrics)
accuracy_nolas <- Accuracy(y_pred_clean, y_true_clean)
f1_nolas <- F1_Score(y_pred_clean, y_true_clean, positive = "Deceased")
balanced_accuracy_nolas <- (Sensitivity(y_pred_clean, y_true_clean, positive = "Deceased") +
                              Specificity(y_pred_clean, y_true_clean, positive = "Deceased")) / 2




cat("📊 Accuracy:", round(accuracy_nolas, 3), "\n")
cat("⚖️ Balanced Accuracy:", round(balanced_accuracy_nolas, 3), "\n")
cat("🎯 F1-Score:", round(f1_nolas, 3), "\n")
cat("📈 AUC:", round(auc_val_nolas, 3), "\n")

# 🧮 Step 26: Compute performance metrics (Accuracy, Balanced Accuracy, F1)
library(MLmetrics)

accuracy <- Accuracy(y_pred_clean, y_true_clean)
f1 <- F1_Score(y_pred_clean, y_true_clean, positive = "Deceased")

bal_acc <- (
  Sensitivity(y_pred_clean, y_true_clean, positive = "Deceased") +
    Specificity(y_pred_clean, y_true_clean, positive = "Deceased")
) / 2

results_nolas <- data.frame(
  Accuracy = accuracy,
  BalancedAccuracy = bal_acc,
  F1 = f1,
  AUC = auc_val_nolas
)
print(results_nolas)

# 🧬 Step 27: Extract most important features (top genes)
# 💬 Description:
# Extract loadings from the NOLAS model for both latent variables (LV1, LV2).
# Select top 150 genes with highest absolute loading for each component.
# 🧬 Step 1: Extract gene loadings from NOLAS result
gene_loadings <- res_nolas_final@loadings

# 🧬 Step 2: Filter RNA features
rna_loadings <- gene_loadings[grepl("^RNASeq", rownames(gene_loadings)), ]
rownames(rna_loadings) <- gsub("^RNASeq\\s+", "", rownames(rna_loadings))

# 🧬 Step 3: Filter miRNA features
mirna_loadings <- gene_loadings[grepl("^miRNA", rownames(gene_loadings)), ]
rownames(mirna_loadings) <- gsub("^miRNA\\s+", "", rownames(mirna_loadings))

# 🧬 Step 4: Select top features from each latent variable for both layers
top_n <- 150  # number of top features per LV

# RNA LV1 + LV2
lv1_rna <- abs(rna_loadings[, "LV1"])
top_rna_lv1 <- names(sort(lv1_rna, decreasing = TRUE))[1:top_n]

lv2_rna <- abs(rna_loadings[, "LV2"])
top_rna_lv2 <- names(sort(lv2_rna, decreasing = TRUE))[1:top_n]

# miRNA LV1 + LV2
lv1_mirna <- abs(mirna_loadings[, "LV1"])
top_mirna_lv1 <- names(sort(lv1_mirna, decreasing = TRUE))[1:top_n]

lv2_mirna <- abs(mirna_loadings[, "LV2"])
top_mirna_lv2 <- names(sort(lv2_mirna, decreasing = TRUE))[1:top_n]

# 🧬 Step 5: Combine and deduplicate
top_genes_all <- unique(c(top_rna_lv1, top_rna_lv2, top_mirna_lv1, top_mirna_lv2))

# 🧮 Final count
length(top_genes_all)
# 🔬 Step 28: Perform functional enrichment analysis using Enrichr
# 💬 Description:
# Analyze the top genes for biological relevance using GO and KEGG pathway databases.
# 📦 Install and load enrichR if not already installed
if (!requireNamespace("enrichR", quietly = TRUE)) install.packages("enrichR")
library(enrichR)

# 🧬 Define databases for enrichment analysis
dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human")

# 🔍 Run enrichment analysis on top genes from RNA + miRNA (NOLAS)
enriched_nolas <- enrichr(top_genes_all, dbs)

# 💾 Save selected genes to a text file
write.table(top_genes_all, "top_genes_nolas_mirna.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 📊 View top enriched GO biological processes
head(enriched_nolas[["GO_Biological_Process_2021"]][, c("Term", "Adjusted.P.value")], 10)

# 📊 View top enriched KEGG pathways
head(enriched_nolas[["KEGG_2021_Human"]][, c("Term", "Adjusted.P.value")], 10)

# 📊 Step 29: Final summary metrics printout
# 💬 Description:
# Print final performance of NOLAS model using standard metrics.
accuracy_nolas <- Accuracy(y_pred_clean, y_true_clean)
f1_nolas <- F1_Score(y_pred_clean, y_true_clean, positive = "Deceased")

balanced_accuracy_nolas <- (
  Sensitivity(y_pred_clean, y_true_clean, positive = "Deceased") +
    Specificity(y_pred_clean, y_true_clean, positive = "Deceased")
) / 2

cat("📊 NOLAS Accuracy:", round(accuracy_nolas, 3), "\n")
cat("⚖️ NOLAS Balanced Accuracy:", round(balanced_accuracy_nolas, 3), "\n")
cat("🎯 NOLAS F1-Score (Deceased):", round(f1_nolas, 3), "\n")

results_nolas <- data.frame(
  Accuracy = accuracy_nolas,
  BalancedAccuracy = balanced_accuracy_nolas,
  F1 = f1_nolas,
  AUC = auc_val_nolas
)
print(results_nolas)

# 🧪 PR-AUC  NOLAS
library(PRROC)
pr_auc_nolas <- pr_curve_nolas$auc.integral
pr_curve_nolas <- pr.curve(
  scores.class0 = y_score_nolas[y_true_clean == "Deceased"],
  scores.class1 = y_score_nolas[y_true_clean == "Alive"],
  curve = TRUE
)

cat("🔬 PR-AUC (NOLAS):", round(pr_curve_nolas$auc.integral, 4), "\n")
#لهون شغال 
# 📊 Score Plot: Latent Variable 1 vs Latent Variable 

# Step 1: Extract scores
scores_nolas <- scores(res_nolas_final)

# Step 2: Match samples with training labels
common_samples <- rownames(scores_nolas)
labels_matched <- labels_train[common_samples]

# Step 3: Prepare data frame for plotting
df_scores <- data.frame(
  LV1 = scores_nolas[, "LV1"],
  LV2 = scores_nolas[, "LV2"],
  Class = factor(labels_matched, levels = c("Alive", "Deceased"))
)

# Step 4: Plot
ggplot(df_scores, aes(x = LV1, y = LV2, color = Class)) +
  geom_point(shape = 16, size = 2) +
  scale_color_manual(values = c("Alive" = "green", "Deceased" = "red")) +
  labs(
    title = "NOLAS Latent Space (LV1 vs LV2)",
    x = "Latent Variable 1 (LV1)",
    y = "Latent Variable 2 (LV2)",
    color = "Survival Status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "right"
  )


# 📈 ROC Curve Plot for NOLAS
# Step 1: Build ROC object
roc_obj_nolas <- roc(
  response = y_true_clean,
  predictor = y_score_nolas,
  levels = c("Alive", "Deceased"),
  direction = "<"
)

# Step 2: Plot ROC curve
plot(roc_obj_nolas, col = "#1c61b6", lwd = 3, main = "ROC Curve - NOLAS Model")
abline(a = 0, b = 1, lty = 2, col = "gray")  # Add diagonal line (random classifier)

# Step 3: Print AUC value
auc_val_nolas <- auc(roc_obj_nolas)
cat("📈 AUC (NOLAS):", round(auc_val_nolas, 4), "\n")
plot(pr_curve_nolas, main = "NOLAS Precision-Recall Curve", color = "blue")
########
# 📦 Load required libraries
library(nolas)
library(caret)

# 🔁 Step 1: Set up Stability Analysis
set.seed(123)  # For reproducibility
n_iterations <- 50  # Number of repetitions
subsample_ratio <- 0.8  # 80% training each time

# 📝 Initialize a list to store selected features for each iteration
selected_features_list <- list()

# 🚀 Step 2: Start the iterations
for (i in 1:n_iterations) {
  
  cat("🔄 Running iteration:", i, "\n")
  
  # ✅ Step 2.1: Create a random 80% subsample of training samples
  train_sub_idx <- createDataPartition(training_metadata$Survival, p = subsample_ratio, list = FALSE)
  
  subsample_train_data <- list(
    RNASeq = omics_train$RNASeq[train_sub_idx, ],
    RPPA   = omics_train$RPPA[train_sub_idx, ],
    miRNA  = omics_train$miRNA[train_sub_idx, ]
  )
  
  subsample_labels <- training_metadata$Survival[train_sub_idx]
  subsample_metadata <- data.frame(Survival = subsample_labels)
  
  # ✅ Step 2.2: Train a NOLAS model on the subsample
  res_nolas_sub <- nolas(
    X = subsample_train_data,
    Y = subsample_metadata,
    nlv = 2,
    imputeNA = TRUE,
    center = TRUE,
    scale = FALSE,
    response = "Survival",
    rvc = 0.25
  )
  
  # ✅ Step 2.3: Extract selected features (non-zero loadings)
  loadings_matrix <- res_nolas_sub@loadings
  selected_features <- rownames(loadings_matrix)[apply(loadings_matrix, 1, function(x) any(x != 0))]
  
  # ✅ Step 2.4: Save selected features for this iteration
  selected_features_list[[i]] <- selected_features
}

# 🧮 Step 3: Summarize stability results
# Flatten the list into a single vector
all_selected_features <- unlist(selected_features_list)

# Count how many times each feature was selected
feature_frequency <- table(all_selected_features)

# Convert to a data frame
stability_results <- data.frame(
  Feature = names(feature_frequency),
  Frequency = as.integer(feature_frequency)
)

# Sort by frequency (descending)
stability_results <- stability_results[order(-stability_results$Frequency), ]

# 📊 Step 4: View top stable features
head(stability_results, 10)

# 💾 Step 5: Save stability results to file
write.csv(stability_results, "stability_results_nolas.csv", row.names = FALSE)

cat("✅ Stability analysis completed and results saved to 'stability_results_nolas.csv'\n")



































#2#############################################
#DIABLO
# ========================
# 1. Prepare Training Data (Filtering)
# ========================
set.seed(123)
# RNA block
rna_df <- as.data.frame(omics_train$RNASeq)
nzv_info_rna <- caret::nearZeroVar(rna_df, saveMetrics = TRUE)
genes_to_keep_rna <- setdiff(colnames(rna_df), rownames(nzv_info_rna[nzv_info_rna$nzv, ]))

# miRNA block
mirna_df <- as.data.frame(omics_train$miRNA)
nzv_info_mirna <- caret::nearZeroVar(mirna_df, saveMetrics = TRUE)
genes_to_keep_mirna <- setdiff(colnames(mirna_df), rownames(nzv_info_mirna[nzv_info_mirna$nzv, ]))

# Build filtered training blocks
omics_train_diablo <- list(
  RNASeq = omics_train$RNASeq[, genes_to_keep_rna, drop = FALSE],
  RPPA   = omics_train$RPPA,
  miRNA  = omics_train$miRNA[, genes_to_keep_mirna, drop = FALSE]
)
omics_test <- list(
  RNASeq = t(multi_omics_test$RNASeq),
  RPPA   = t(multi_omics_test$RPPA),
  miRNA  = t(multi_omics_test$miRNA)
)

# ========================
# 2. Prepare Testing Data (Filtering)
# ========================

omics_test$RNASeq <- t(multi_omics_test$RNASeq)
omics_test$miRNA  <- t(multi_omics_test$miRNA)
omics_test$RPPA   <- t(multi_omics_test$RPPA)

common_genes_rna   <- intersect(genes_to_keep_rna, colnames(omics_test$RNASeq))
common_genes_mirna <- intersect(genes_to_keep_mirna, colnames(omics_test$miRNA))

omics_test_diablo <- list(
  RNASeq = omics_test$RNASeq[, common_genes_rna, drop = FALSE],
  RPPA   = omics_test$RPPA,
  miRNA  = omics_test$miRNA[, common_genes_mirna, drop = FALSE]
)

# Align samples
common_test_samples <- Reduce(intersect, list(
  rownames(omics_test_diablo$RNASeq),
  rownames(omics_test_diablo$RPPA),
  rownames(omics_test_diablo$miRNA)
))

omics_test_diablo <- lapply(omics_test_diablo, function(block) block[common_test_samples, , drop = FALSE])

# ========================
# 3. Train DIABLO Model
# ========================

survival_labels_train <- labels[rownames(omics_train_diablo$RNASeq)]

keepX_list <- list(
  RNASeq = rep(80, 2),
  RPPA   = rep(20, 2),
  miRNA  = rep(20, 2)
)

design_matrix <- matrix(0.6, nrow = 3, ncol = 3)
diag(design_matrix) <- 0
rownames(design_matrix) <- colnames(design_matrix) <- names(omics_train_diablo)

library(mixOmics)
diablo_model <- block.splsda(
  X = omics_train_diablo,
  Y = survival_labels_train,
  ncomp = 2,
  design = design_matrix,
  keepX = keepX_list
)

# ========================
# 4. Predict on Test Set
# ========================

predicted_labels <- pred_diablo$MajorityVote$max.dist[, "comp2"]


common_ids <- intersect(names(predicted_labels), names(labels_test))
predicted_labels <- predicted_labels[common_ids]
true_labels <- labels_test[common_ids]


predicted_labels <- factor(predicted_labels, levels = c("Alive", "Deceased"))
true_labels <- factor(true_labels, levels = c("Alive", "Deceased"))


conf_mat <- caret::confusionMatrix(predicted_labels, true_labels)
print(conf_mat)

# 5. Evaluate Performance
# ========================

conf_mat <- caret::confusionMatrix(predicted_labels, true_labels)
conf_mat <- caret::confusionMatrix(predicted_labels, true_labels)
conf_mat
accuracy_diablo          <- conf_mat$overall["Accuracy"]
balanced_accuracy_diablo <- conf_mat$byClass["Balanced Accuracy"]
f1_diablo                <- conf_mat$byClass["F1"]

# ROC AUC for each block
roc_rna   <- roc(true_labels, pred_diablo$predict$RNASeq[, , 2][, "Deceased"], levels = c("Alive", "Deceased"))
roc_mirna <- roc(true_labels, pred_diablo$predict$miRNA[, , 2][, "Deceased"], levels = c("Alive", "Deceased"))
roc_rppa  <- roc(true_labels, pred_diablo$predict$RPPA[, , 2][, "Deceased"], levels = c("Alive", "Deceased"))

auc_rna   <- auc(roc_rna)
auc_mirna <- auc(roc_mirna)
auc_rppa  <- auc(roc_rppa)
auc_mean  <- mean(c(auc_rna, auc_mirna, auc_rppa))

# PR-AUC
pr_curve_diablo <- pr.curve(
  scores.class0 = pred_diablo$predict$RNASeq[, , 2][true_labels == "Deceased", "Deceased"],
  scores.class1 = pred_diablo$predict$RNASeq[, , 2][true_labels == "Alive", "Deceased"],
  curve = TRUE
)

pr_auc_diablo <- pr_curve_diablo$auc.integral

# ========================
# 6. Plot Score Plot (LV1 vs LV2)
# ========================

survival_labels_train_fixed <- factor(survival_labels_train[rownames(diablo_model$X$RNASeq)], levels = c("Alive", "Deceased"))

plotIndiv(
  object = diablo_model,
  comp = c(1, 2),
  group = survival_labels_train_fixed,
  ind.names = FALSE,
  legend = TRUE,
  title = "DIABLO Score Plot (Comp1 vs Comp2)",
  col = c("green3", "red3")
)

# ========================
# 7. Stability Analysis (50x Subsampling)
# ========================

library(mixOmics)
library(caret)

# 🧪 Parameters
n_iterations <- 50
set.seed(123)  # لتكرار النتائج

# 🧪 Initialize list to collect selected genes
selected_genes_all <- list()

for (i in 1:n_iterations) {
  cat("Iteration:", i, "\n")
  
  # 🧪 Subsample 80% of training data
  set.seed(i)  # Different seed per iteration
  sample_idx <- createDataPartition(survival_labels_train, p = 0.8, list = FALSE)
  
  # 🧪 Subset training data
  omics_train_subsample <- lapply(omics_train_diablo, function(x) x[sample_idx, , drop = FALSE])
  labels_subsample <- survival_labels_train[sample_idx]
  
  # 🧪 Train DIABLO model on subsample
  diablo_subsample_model <- block.splsda(
    X = omics_train_subsample,
    Y = labels_subsample,
    ncomp = 2,
    design = design_matrix,
    keepX = keepX_list
  )
  
  # 🧪 Extract selected features
  genes_rna_lv1   <- selectVar(diablo_subsample_model, comp = 1)$RNASeq$name
  genes_rna_lv2   <- selectVar(diablo_subsample_model, comp = 2)$RNASeq$name
  mirna_lv1       <- selectVar(diablo_subsample_model, comp = 1)$miRNA$name
  mirna_lv2       <- selectVar(diablo_subsample_model, comp = 2)$miRNA$name
  
  # 🧪 Combine all selected genes
  selected_genes <- unique(c(genes_rna_lv1, genes_rna_lv2, mirna_lv1, mirna_lv2))
  
  selected_genes_all[[i]] <- selected_genes
}

# 🧮 Step 2: Aggregate selected genes
all_selected_genes <- unlist(selected_genes_all)

# 🧮 Step 3: Count frequency of each gene
gene_frequency <- sort(table(all_selected_genes), decreasing = TRUE)

# 📊 View top selected genes
head(gene_frequency, 10)

# 🧬 Select genes appearing in at least 40 iterations
stable_genes_diablo <- names(gene_frequency[gene_frequency >= 40])

cat("🔬 Number of stable genes (appearing ≥40 times):", length(stable_genes_diablo), "\n")

# 📄 Save stable genes if needed
write.table(stable_genes_diablo, "stable_genes_diablo.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
# ========================
# 8. Functional Enrichment Analysis
# ========================
# 📦 Install enrichR if not installed
if (!requireNamespace("enrichR", quietly = TRUE)) install.packages("enrichR")
library(enrichR)

# 🧬 Databases to use
dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human")

# 🧬 Perform enrichment
enriched_diablo <- enrichr(stable_genes_diablo, dbs)

# 📊 Show top enriched GO terms
head(enriched_diablo[["GO_Biological_Process_2021"]][, c("Term", "Adjusted.P.value")], 10)

# 📊 Show top enriched KEGG pathways
head(enriched_diablo[["KEGG_2021_Human"]][, c("Term", "Adjusted.P.value")], 10)

# ========================
# 9. Final Metrex
metrics_diablo <- list(
  Accuracy          = accuracy_diablo,
  BalancedAccuracy  = balanced_accuracy_diablo,
  F1                = f1_diablo,
  AUC_RNA           = auc_rna,
  AUC_miRNA         = auc_mirna,
  AUC_RPPA          = auc_rppa,
  AUC_mean          = auc_mean,
  PR_AUC            = pr_auc_diablo
)

print(metrics_diablo)




# 📌 Step 1: Prepare matching sample IDs
common_ids <- Reduce(intersect, list(names(predicted_labels), names(y_pred_clean), names(y_true_clean)))

# 📌 Step 2: Align predictions and true labels
pred_diablo_aligned <- predicted_labels[common_ids]  # DIABLO predictions
pred_nolas_aligned  <- y_pred_clean[common_ids]       # NOLAS predictions
truth_aligned       <- y_true_clean[common_ids]       # True labels

# 📌 Step 3: Convert predictions to correct/incorrect (logical vectors)
correct_diablo <- pred_diablo_aligned == truth_aligned
correct_nolas  <- pred_nolas_aligned == truth_aligned

# 📌 Step 4: Build the contingency table
tbl <- table(DIABLO = correct_diablo, NOLAS = correct_nolas)

# 📌 Step 5: Apply McNemar's Test
mcnemar_result <- mcnemar.test(tbl)
print(mcnemar_result)



# 📌 Step: Create a summary performance table comparing DIABLO and NOLAS

# ✅ Build performance table
performance_summary <- data.frame(
  Model = c("DIABLO", "NOLAS"),
  Accuracy = c(metrics_diablo$Accuracy, results_nolas$Accuracy),
  BalancedAccuracy = c(metrics_diablo$BalancedAccuracy, results_nolas$BalancedAccuracy),
  F1 = c(metrics_diablo$F1, results_nolas$F1),
  AUC = c(metrics_diablo$AUC_mean, results_nolas$AUC),
  PR_AUC = c(metrics_diablo$PR_AUC, pr_curve_nolas$auc.integral)
)

# ✅ Print the final table
print(performance_summary)

