# =================================================================
# COMPREHENSIVE CMML PREDICTIVE MODELING PIPELINE
# =================================================================

# PURPOSE: This script implements a comprehensive analysis pipeline for
# CMML (Chronic Myelomonocytic Leukemia) classification, including:
# 1. Data processing and feature engineering
# 2. Unsupervised clustering with multiple methods
# 3. Cluster visualization with t-SNE
# 4. Multiple supervised learning models
# 5. Penalized regression (LASSO, Elastic Net, MCP) with inference
# 6. Model evaluation and comparison
# 7. Prediction for OM-CMML patients


# TODO: The covariates are all indicators, we should treat them as binary variable instead of continuous;
# In elastic net, I think we should not standardize the data
# For random forest, add SHARP values, which makes more sense


# ----------------------------------------------------------------------
# SECTION 0: Loading Required Packages
# ----------------------------------------------------------------------

cat("\n==== CHECKING AND LOADING REQUIRED PACKAGES ====\n")

# Core packages for data manipulation, clustering, and visualization
core_packages <- c(
  "tidyverse",     # Collection of packages for data manipulation and visualization
  "vegan",         # Ordination analysis and distance calculations
  "cluster",       # Clustering methods 
  "mclust",        # Model-based clustering
  "Rtsne",         # t-SNE visualization
  "cowplot",       # Enhanced ggplot layouts
  "dbscan",        # Density-based clustering
  "kernlab",       # Kernel-based machine learning
  "ggplot2",       # Plotting library
  "reshape2",      # Data reshaping for visualization
  "gridExtra"      # Arrange multiple plots on a page
)

# Packages for supervised learning
ml_packages <- c(
  "caret",         # Framework for training machine learning models
  "randomForest",  # Implementation of Random Forest algorithm
  "e1071",         # Contains support vector machines implementation
  "nnet",          # Neural network capabilities
  "xgboost",       # Extreme Gradient Boosting algorithm
  "glmnet"         # Implements elastic net, ridge, and lasso regression
)

# Packages for statistical inference
inference_packages <- c(
  "ncvreg",          # MCP regression
  "hdm",             # High-dimensional metrics package for debiased Lasso
  "hdi",             # High-dimensional inference
  "selectiveInference", # Post-selection inference for Lasso
  "boot",            # Bootstrap methods
  "PRROC",           # Precision-recall and ROC curves
  "pROC"             # Tools for ROC curve analysis
)

# Combine all required packages
all_packages <- c(core_packages, ml_packages, inference_packages)

# This loop checks each package and installs it if it's not already installed
for(pkg in all_packages){
  if (!require(pkg, character.only = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  } else {
    cat("Package already installed:", pkg, "\n")
  }
}

# ----------------------------------------------------------------------
# SECTION 1: Loading and Preparing the Data
# ----------------------------------------------------------------------

cat("\n==== LOADING AND PREPARING DATA ====\n")

# Set the working directory - MODIFY THIS PATH TO MATCH YOUR SYSTEM
# setwd("/your/path/here")  # Uncomment and modify this line

# Load the data file
cat("Loading data file...\n")
data <- read.csv("OM-CMML_for_analysis_April17.csv", header = TRUE, stringsAsFactors = FALSE)

# Display data structure to verify successful loading
cat("Examining data structure:\n")
str(data)

# Summarize the data
cat("\nData summary:\n")
cat("Total number of patients:", nrow(data), "\n")
cat("Number of variables:", ncol(data), "\n")
diagnostic_counts <- table(data$Diagnostic.group)
cat("Patients by diagnostic group:\n")
print(diagnostic_counts)

# ------------------------------
# 1.1 Process and Convert Selected Columns 
# ------------------------------
# This section extracts and processes the columns needed for prediction

cat("\nProcessing features for analysis...\n")

# Define columns for clustering/predictors
# Assuming specific column structure (modify as needed for your dataset)
selected_cols <- c(7:ncol(data))


# Process clinical columns (typically columns 4-6)
# Convert Yes/No responses to numeric (1/0)
convert_yes_no <- function(x) {
  x <- tolower(as.character(x))
  as.numeric(ifelse(x == "yes", 1, ifelse(x == "no", 0, NA)))
}


# Process binary columns (typically columns 7 to end)
binary_cols <- data[, 7:ncol(data)]
binary_converted <- as.data.frame(lapply(binary_cols, function(x) {
  x <- tolower(as.character(x))
  as.numeric(ifelse(x == "yes", 1, ifelse(x == "no", 0, NA)))
}))

# ------------------------------
# 1.2 Handle Missing Values 
# ------------------------------

cat("\nImputing missing values...\n")

# Function to impute missing values using mode
impute_mode <- function(x) {
  # Ensure x is numeric
  x <- as.numeric(x)
  
  # Only proceed if there are missing values
  if (any(is.na(x))) {
    # Create a frequency table of non-NA values
    freq <- table(x)
    
    # Identify the mode(s) - values with the maximum frequency
    mode_vals <- as.numeric(names(freq)[freq == max(freq)])
    
    # In case of ties, choose the first one
    mode_val <- mode_vals[1]
    
    # Replace NA values with the mode
    x[is.na(x)] <- mode_val
  }
  return(x)
}

# Apply imputation to all numeric columns
binary_converted <- as.data.frame(lapply(binary_converted, impute_mode))

bmmono = data$BM_Mono
bmmono[is.na(bmmono)] = mean(bmmono, na.rm = TRUE)
bmmono = scale(bmmono)
bmmono = as.data.frame(bmmono)
colnames(bmmono) = c("standardized_BM_Mono")
# ------------------------------
# 1.3 Combine All Features 
# ------------------------------

# Combine all processed features into a single matrix
clust_data <- bind_cols(binary_converted,bmmono)

# Display feature matrix structure
cat("\nFeature matrix structure:\n")
str(clust_data)

# ----------------------------------------------------------------------
# SECTION 2: Unsupervised Clustering Methods
# ----------------------------------------------------------------------

cat("\n==== UNSUPERVISED CLUSTERING METHODS ====\n")

# Create a list to store all clustering results
cluster_results <- list()


# ------------------------------
# 2.2 K-means Clustering
# ------------------------------
cat("\nPerforming K-means clustering...\n")

# Set seed for reproducibility
set.seed(123)

# Perform K-means clustering with k=10
kmeans_res <- kmeans(clust_data, centers = 10, nstart = 25)

# Store clustering results
data$kmeans_cluster <- as.factor(kmeans_res$cluster)
cluster_results$kmeans <- kmeans_res$cluster

# Display cross-tabulation
cat("K-means Clustering Cross-tabulation:\n")
print(table(DiagnosticGroup = data$Diagnostic.group, KMeans_Cluster = data$kmeans_cluster))

# ------------------------------
# 2.3 PAM Clustering
# ------------------------------
cat("\nPerforming PAM (Partitioning Around Medoids) clustering...\n")

# Perform PAM clustering with k=3
pam_res <- pam(clust_data, k = 10)

# Store clustering results
data$pam_cluster <- as.factor(pam_res$clustering)
cluster_results$pam <- pam_res$clustering

# Display cross-tabulation
cat("PAM Clustering Cross-tabulation:\n")
print(table(DiagnosticGroup = data$Diagnostic.group, PAM_Cluster = data$pam_cluster))

# ------------------------------
# 2.4 Model-based Clustering (Mclust)
# ------------------------------
cat("\nPerforming model-based clustering...\n")

# Perform model-based clustering
mclust_res <- Mclust(clust_data)

# Store clustering results
data$mclust_cluster <- as.factor(mclust_res$classification)
cluster_results$mclust <- mclust_res$classification

# Display cross-tabulation
cat("Mclust Clustering Cross-tabulation:\n")
print(table(DiagnosticGroup = data$Diagnostic.group, Mclust_Cluster = data$mclust_cluster))

# ------------------------------
# 2.5 Spectral Clustering
# ------------------------------
cat("\nPerforming spectral clustering...\n")

# Perform spectral clustering
spec_res <- specc(as.matrix(clust_data), centers = 10)

# Store clustering results
data$spectral_cluster <- as.factor(spec_res)
cluster_results$spectral <- spec_res

# Display cross-tabulation
cat("Spectral Clustering Cross-tabulation:\n")
print(table(DiagnosticGroup = data$Diagnostic.group, Spectral_Cluster = data$spectral_cluster))

# ------------------------------
# 2.6 Save All Clustering Results
# ------------------------------
cat("\nSaving clustering results...\n")

# Create a data frame with all clustering results
clustering_df <- data.frame(
  PatientID = data$Patient.ID,
  DiagnosticGroup = data$Diagnostic.group,
  KMeans_Cluster = data$kmeans_cluster,
  PAM_Cluster = data$pam_cluster,
  Mclust_Cluster = data$mclust_cluster,
  Spectral_Cluster = data$spectral_cluster
)

# Save to CSV
write.csv(clustering_df, "CMML_Clustering_Results.csv", row.names = FALSE)

# ----------------------------------------------------------------------
# SECTION 3: Visualization with t-SNE
# ----------------------------------------------------------------------

cat("\n==== VISUALIZING CLUSTERS WITH T-SNE ====\n")

# ------------------------------
# 3.1 Run t-SNE
# ------------------------------
cat("\nRunning t-SNE dimensionality reduction...\n")

# Add small jitter to avoid duplicates
set.seed(123)
clust_data_jitter <- clust_data + matrix(rnorm(nrow(clust_data) * ncol(clust_data), 
                                         mean = 0, sd = 1e-8),
                                         nrow = nrow(clust_data))

# Run t-SNE
tsne_out <- Rtsne(as.matrix(clust_data_jitter), dims = 2, perplexity = 30,
                  verbose = TRUE, max_iter = 500, check_duplicates = FALSE)

# Create t-SNE data frame for plotting
tsne_df <- data.frame(
  tSNE1 = tsne_out$Y[,1],
  tSNE2 = tsne_out$Y[,2],
  DiagnosticGroup = data$Diagnostic.group,
  KMeans_Cluster = data$kmeans_cluster,
  PAM_Cluster = data$pam_cluster,
  Mclust_Cluster = data$mclust_cluster,
  Spectral_Cluster = data$spectral_cluster
)

# ------------------------------
# 3.2 Create t-SNE Visualizations
# ------------------------------
cat("\nCreating t-SNE visualizations for each clustering method...\n")

# Define publication-quality theme
pub_theme <- theme_classic(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5))

# Create t-SNE plots for each clustering method

# K-means clustering
p2 <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = KMeans_Cluster, shape = DiagnosticGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE: K-means Clusters", x = "t-SNE Dim 1", y = "t-SNE Dim 2") +
  pub_theme

# PAM clustering
p3 <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = PAM_Cluster, shape = DiagnosticGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE: PAM Clusters", x = "t-SNE Dim 1", y = "t-SNE Dim 2") +
  pub_theme

# Model-based clustering
p4 <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Mclust_Cluster, shape = DiagnosticGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE: Mclust Clusters", x = "t-SNE Dim 1", y = "t-SNE Dim 2") +
  pub_theme

# Spectral clustering
p5 <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Spectral_Cluster, shape = DiagnosticGroup)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE: Spectral Clusters", x = "t-SNE Dim 1", y = "t-SNE Dim 2") +
  pub_theme



pdf("CMML_tSNE_KMeans_Clusters.pdf", width = 8, height = 7)
print(p2)
dev.off()

pdf("CMML_tSNE_PAM_Clusters.pdf", width = 8, height = 7)
print(p3)
dev.off()

pdf("CMML_tSNE_Mclust_Clusters.pdf", width = 8, height = 7)
print(p4)
dev.off()

pdf("CMML_tSNE_Spectral_Clusters.pdf", width = 8, height = 7)
print(p5)
dev.off()

# Combine all plots into a single figure
pdf("CMML_tSNE_All_Clusters.pdf", width = 16, height = 12)
grid.arrange(p2, p3, p4, p5, ncol = 2)
dev.off()

cat("t-SNE visualizations saved successfully.\n")

# ----------------------------------------------------------------------
# SECTION 4: Preparing for Supervised Learning
# ----------------------------------------------------------------------
cat("\n==== SETTING UP SUPERVISED LEARNING ====\n")

# Define training set (MD-CMML and MDS) and test set (OM-CMML)
train_index <- which(data$Diagnostic.group %in% c("MD-CMML", "MDS"))
test_index <- which(data$Diagnostic.group == "OM-CMML")

# Extract features for training and testing
train_X <- clust_data[train_index, ]
test_X <- clust_data[test_index, ]

# Ensure we have binary classification (MD-CMML vs MDS)
# Convert response to a two-level factor with VALID R variable names
# Replace hyphens with underscores to make valid R variable names
train_Y_orig <- data$Diagnostic.group[train_index]
train_Y <- factor(ifelse(train_Y_orig == "MD-CMML", "MD_CMML", "MDS"))

# Define positive class for consistent evaluation
positive_class <- "MD_CMML"

# Display class distribution in training data
cat("\nTraining data class distribution:\n")
print(table(train_Y))
cat("Number of training samples:", length(train_Y), "\n")
cat("Number of features:", ncol(train_X), "\n")

# Set up cross-validation for binary classification
set.seed(123)
train_control <- trainControl(
  method = "cv",             # Use cross-validation
  number = 10,               # 10-fold cross-validation
  search = "grid",           # Use grid search for parameter tuning
  savePredictions = "final", # Save predictions from cross-validation
  classProbs = TRUE,         # Calculate class probabilities (required for ROC)
  summaryFunction = twoClassSummary, # Use metrics for binary classification
  verboseIter = TRUE,        # Display progress during training
  allowParallel = TRUE       # Allow parallel processing to speed up computations
)

cat("Cross-validation configured with 10 folds for binary classification.\n")

# ----------------------------------------------------------------------
# SECTION 5: Training Classification Models with Inference
# ----------------------------------------------------------------------

cat("\n==== TRAINING CLASSIFICATION MODELS WITH INFERENCE ====\n")

# Create a list to store all trained models
model_list <- list()

# ------------------------------
# 5.2 Elastic Net with Bootstrap Inference
# ------------------------------
cat("\nTraining Elastic Net model with bootstrap inference...\n")
set.seed(123)

# Train Elastic Net model with standard parameter tuning
enet_grid <- expand.grid(
  alpha = 0.5,      # From ridge (0) to lasso (1)
  lambda = 10^seq(-5, 0, length = 100) # Range of regularization parameters
)

enet_model <- train(
  x = train_X,
  y = train_Y,
  method = "glmnet",        # Elastic Net
  trControl = train_control,
  tuneGrid = enet_grid,     # Only specify alpha grid, let lambda be automatic
  metric = "ROC",
  preProcess = NULL,        # Disable preprocessing/standardization
  standardize = FALSE       # Disable standardization within glmnet
)
model_list$ElasticNet <- enet_model

# Display model performance
cat("Elastic Net - Best Tuning Parameters:\n")
print(enet_model$bestTune)
cat("Elastic Net - Cross-Validation Performance (Best Parameters):\n")
print(enet_model$results[enet_model$results$alpha == enet_model$bestTune$alpha, ][1,])
enet_pred <- predict(enet_model, newdata = train_X)
cat("Elastic Net - Training Confusion Matrix:\n")
print(confusionMatrix(enet_pred, train_Y, positive = positive_class))

# Extract coefficients from the final model
enet_coef <- as.vector(coef(enet_model$finalModel, s = enet_model$bestTune$lambda))
enet_coef_names <- rownames(coef(enet_model$finalModel, s = enet_model$bestTune$lambda))
enet_coef_df <- data.frame(
  Variable = enet_coef_names[-1],  # Remove intercept
  Coefficient = enet_coef[-1]      # Remove intercept
)
enet_coef_df <- enet_coef_df[order(abs(enet_coef_df$Coefficient), decreasing = TRUE), ]

# Save non-zero coefficients
nonzero_enet <- enet_coef_df[enet_coef_df$Coefficient != 0, ]
write.csv(nonzero_enet, "CMML_ElasticNet_NonZero_Coefficients.csv", row.names = FALSE)

# Create coefficient plot
pdf("CMML_ElasticNet_Coefficients.pdf", width = 10, height = 8)
ggplot(head(nonzero_enet, 20), aes(x = reorder(Variable, abs(Coefficient)), y = Coefficient)) +
  geom_bar(stat = "identity", aes(fill = Coefficient > 0)) +
  scale_fill_manual(values = c("firebrick", "dodgerblue"), 
                    labels = c("Negative", "Positive"),
                    name = "Direction") +
  coord_flip() +
  labs(title = "Elastic Net - Top 20 Variables by Coefficient Magnitude",
       x = "Variable",
       y = "Coefficient Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
dev.off()

# ------------------------------
# 5.2.1 Non-parametric Bootstrap for Elastic Net
# ------------------------------
cat("\nPerforming non-parametric bootstrap for Elastic Net inference...\n")

# Number of bootstrap iterations
n_bootstrap <- 500
X_matrix <- as.matrix(train_X)
bootstrap_results <- matrix(NA, nrow = n_bootstrap, ncol = ncol(X_matrix))
colnames(bootstrap_results) <- colnames(X_matrix)

set.seed(42)
for(i in 1:n_bootstrap) {
  if(i %% 50 == 0) cat("  Non-parametric bootstrap iteration:", i, "\n")
  
  # Sample with replacement
  boot_idx <- sample(nrow(X_matrix), replace = TRUE)
  X_boot <- X_matrix[boot_idx, ]
  y_boot <- train_Y[boot_idx]
  
  # Fit Elastic Net with fixed parameters from the original model
  tryCatch({
    # Use caret's train for consistency with original model
    boot_model <- train(
      x = X_boot,
      y = y_boot,
      method = "glmnet",
      trControl = trainControl(method = "none"),
      tuneGrid = data.frame(
        alpha = enet_model$bestTune$alpha,
        lambda = enet_model$bestTune$lambda
      ),
      preProcess = NULL,        # Disable preprocessing/standardization
      standardize = FALSE       # Disable standardization within glmnet
    )
    
    # Extract coefficients
    coef_boot <- as.vector(coef(boot_model$finalModel, s = boot_model$bestTune$lambda))[-1]
    bootstrap_results[i, ] <- coef_boot
  }, error = function(e) {
    cat("Warning: Error in bootstrap iteration", i, ":", conditionMessage(e), "\n")
  })
}

# Calculate 95% confidence intervals
ci_lower <- apply(bootstrap_results, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
ci_upper <- apply(bootstrap_results, 2, function(x) quantile(x, 0.975, na.rm = TRUE))

# Calculate bootstrap p-values (proportion of bootstrap samples where coefficient has opposite sign)
p_values <- sapply(1:ncol(bootstrap_results), function(j) {
  coef_val <- enet_coef_df$Coefficient[enet_coef_df$Variable == colnames(bootstrap_results)[j]]
  if (coef_val == 0) return(1) # If coefficient is 0, p-value is 1
  opposite_sign <- if (coef_val > 0) {
    bootstrap_results[, j] < 0
  } else {
    bootstrap_results[, j] > 0
  }
  mean(opposite_sign, na.rm = TRUE) * 2 # Two-sided p-value
})

# Determine significance based on confidence intervals and p-values
ci_significant <- (ci_lower > 0 & ci_upper > 0) | (ci_lower < 0 & ci_upper < 0)
p_significant <- p_values < 0.05

# Create bootstrap results data frame
bootstrap_inference <- data.frame(
  Variable = colnames(X_matrix),
  Coefficient = enet_coef[-1],
  CI_Lower = ci_lower,
  CI_Upper = ci_upper,
  P_Value = p_values,
  CI_Significant = ci_significant,
  P_Significant = p_significant
)

# Sort by absolute coefficient value
bootstrap_inference <- bootstrap_inference[order(abs(bootstrap_inference$Coefficient), decreasing = TRUE), ]

# Save significant variables
significant_enet_vars <- bootstrap_inference[bootstrap_inference$CI_Significant, ]
write.csv(significant_enet_vars, "CMML_ElasticNet_Significant_Variables.csv", row.names = FALSE)

# Create plot of significant variables with confidence intervals
pdf("CMML_ElasticNet_Bootstrap_Inference.pdf", width = 10, height = 8)
ggplot(head(significant_enet_vars, 20), 
       aes(x = reorder(Variable, abs(Coefficient)), y = Coefficient)) +
  geom_bar(stat = "identity", aes(fill = Coefficient > 0)) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  scale_fill_manual(values = c("firebrick", "dodgerblue"), 
                    labels = c("Negative", "Positive"),
                    name = "Direction") +
  coord_flip() +
  labs(title = "Elastic Net - Significant Variables with Bootstrap Confidence Intervals",
       x = "Variable",
       y = "Coefficient Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
dev.off()


# ------------------------------
# 5.4 Generate and Save Predictions for Testing Data
# ------------------------------
cat("\nGenerating predictions for OM-CMML test data...\n")

# Create a data frame to store predictions
test_predictions <- data.frame(
  PatientID = data$Patient.ID[test_index],
  ActualGroup = data$Diagnostic.group[test_index]
)

# Add Elastic Net predictions
enet_test_pred <- predict(enet_model, newdata = test_X)
enet_test_prob <- predict(enet_model, newdata = test_X, type = "prob")
test_predictions$ElasticNet_Prediction <- enet_test_pred
# Convert prediction back to original class name format
test_predictions$ElasticNet_Prediction_Original <- ifelse(test_predictions$ElasticNet_Prediction == positive_class, 
                                                         "MD-CMML", "MDS")
test_predictions$ElasticNet_Probability_MD_CMML <- enet_test_prob[, positive_class]



# Add risk score and category based on average probability
test_predictions$Risk_Score <- test_predictions$ElasticNet_Probability_MD_CMML

test_predictions$Risk_Category <- cut(
  test_predictions$Risk_Score,
  breaks = c(0, 0.25, 0.5, 0.75, 1),
  labels = c("Low", "Moderate", "High", "Very High"),
  include.lowest = TRUE
)

# Save predictions to CSV
write.csv(test_predictions, "CMML_Test_Predictions.csv", row.names = FALSE)
