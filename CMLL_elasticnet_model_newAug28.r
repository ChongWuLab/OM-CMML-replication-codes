
library(tidyverse)   # For data manipulation and visualization
library(glmnet)      # For elastic net regression
library(caret)       # For model training and evaluation
library(pROC)        # For ROC curve analysis
library(OptimalCutpoints) # For finding optimal cutoffs

# Set seed for reproducibility
set.seed(123)

# -----------------------------------------------
# 1. Load and Examine Data
# -----------------------------------------------
cat("\n==== LOADING AND EXAMINING DATA ====\n")

# Load the data
data <- read.csv("OM-CMML_for_analysis_for_model.csv", stringsAsFactors = FALSE)

# Basic data examination
cat("Total number of patients:", nrow(data), "\n")
cat("Number of variables:", ncol(data), "\n")

# Check the outcome variable
cat("\nCMML progression distribution:\n")
print(table(data$CMML_progression))

# -----------------------------------------------
# 2. Data Preprocessing
# -----------------------------------------------
cat("\n==== PREPROCESSING DATA ====\n")

# Convert "yes"/"no" to 1/0 for binary variables including the outcome
convert_binary <- function(x) {
  x <- tolower(as.character(x))
  as.numeric(ifelse(x == "yes", 1, ifelse(x == "no", 0, NA)))
}

# Identify binary columns (assuming all except ID columns and BM_Mono are binary)
id_cols <- c("Institution_ID", "Number_ID", "Patient_ID", "Diagnostic_group")
numeric_cols <- c("BM_Mono")
binary_cols <- setdiff(names(data), c(id_cols, numeric_cols))

# Convert binary columns
for (col in binary_cols) {
  data[[col]] <- convert_binary(data[[col]])
}

# Ensure BM_Mono is numeric
data$BM_Mono <- as.numeric(data$BM_Mono)

# Create separate data frames for predictors and outcome
predictors <- data %>% 
  select(-c(id_cols, "CMML_progression")) 

predictors = predictors[,colMeans(predictors,na.rm = TRUE)>0.036]
indx = colSums(is.na(predictors))>30
predictors = predictors[,!indx]

indx = rowSums(is.na(predictors))>10
predictors = predictors[!indx,]

## remove SRSF2 due to high correlation with SRSF2.TET2
predictors = predictors[,!colnames(predictors) %in% "SRSF2"]

outcome <- data$CMML_progression

outcome = outcome[!indx]

# -----------------------------------------------
# 3. Imputation using Mode
# -----------------------------------------------
cat("\n==== PERFORMING MODE IMPUTATION ====\n")

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

# Apply imputation to binary columns
predictors_imputed <- as.data.frame(lapply(predictors, function(x) {
  if(any(is.na(x))) {
    return(impute_mode(x))
  } else {
    return(x)
  }
}))

# For BM_Mono, impute with mean if it exists in the dataset
if("BM_Mono" %in% colnames(predictors_imputed)) {
  bm_mono_col <- predictors_imputed$BM_Mono
  if(any(is.na(bm_mono_col))) {
    predictors_imputed$BM_Mono[is.na(bm_mono_col)] <- mean(bm_mono_col, na.rm = TRUE)
  }
}

# Check if imputation was successful
cat("Number of missing values after imputation:", sum(is.na(predictors_imputed)), "\n")

# -----------------------------------------------
# 4. Find Optimal BM_Mono Cutoff
# -----------------------------------------------
cat("\n==== FINDING OPTIMAL BM_MONO CUTOFF ====\n")

# Check if BM_Mono exists in the dataset
if("BM_Mono" %in% colnames(predictors_imputed)) {
  # Prepare data for cutpoint analysis
  cutpoint_data <- data.frame(
    BM_Mono = predictors_imputed$BM_Mono,
    CMML_progression = factor(outcome, levels = c(0, 1), 
                            labels = c("no", "yes"))
  )
  
  # Find optimal cutpoint using Youden's index
  cut_result <- optimal.cutpoints(
    X = "BM_Mono",
    status = "CMML_progression",
    tag.healthy = "no",
    methods = c("Youden", "SpEqualSe"),
    data = cutpoint_data
  )
  
  # Extract cutpoints
  youden_cutpoint <- cut_result$Youden$Global$optimal.cutoff$cutoff
  speqse_cutpoint <- cut_result$SpEqualSe$Global$optimal.cutoff$cutoff
  
  cat("Optimal BM_Mono cutpoints:\n")
  cat("- Youden index method:", youden_cutpoint, "\n")
  cat("- Sensitivity = Specificity method:", speqse_cutpoint, "\n")
  
  # Create binary variable for BM_Mono based on optimal cutpoint
  predictors_imputed$BM_Mono_Binary <- as.numeric(predictors_imputed$BM_Mono >= youden_cutpoint)
  
  cat("Created binary BM_Mono variable with cutoff:", youden_cutpoint, "\n")
}

predictors_imputed = predictors_imputed[,!colnames(predictors_imputed) %in% c("SRSF2.TET2.ASXL1","TET2","NRAS","KRAS","CBL","BRAF","NF1","PTPN11","BM_Mono")]
# -----------------------------------------------
# 5. Elastic Net Model Training
# -----------------------------------------------
cat("\n==== TRAINING ELASTIC NET MODEL ====\n")

# Prepare data for modeling
X <- as.matrix(predictors_imputed)
y <- outcome


# Set up cross-validation
train_control <- trainControl(
  method = "cv",             # Cross-validation
  number = 10,               # 10-fold cross-validation
  classProbs = TRUE,         # Calculate probabilities
  summaryFunction = twoClassSummary,  # For binary classification
  savePredictions = "final",         # Save final predictions
  verboseIter = TRUE         # Display progress during training
)

# Define parameter grid for elastic net
enet_grid <- expand.grid(
  alpha = seq(0, 1, by = 0.1),     # From ridge (0) to lasso (1)
  lambda = 10^seq(-5, 0, length = 100) # Range of regularization parameters
)


cost_matrix <- matrix(c(0, 10, 1, 0), nrow = 2)

# For glmnet, you can implement this by reweighting observations
# Higher weight for negative cases = higher penalty for misclassifying them
weights <- ifelse(y == 0, cost_matrix[1,2], cost_matrix[2,1])

penalty_factors <- rep(1, ncol(X))  # Initialize all variables with penalty = 1
#penalty_factors[3] <- 0  # Set the variable you don't want penalized to 0



# Train elastic net model
enet_model <- train(
  x = X,
  y = factor(y, levels = c(0, 1), labels = c("no", "yes")),
  method = "glmnet",        # Elastic Net
  trControl = train_control,
  tuneGrid = enet_grid,     
  metric = "Spec",
  family = "binomial",      # For binary classification
  standardize = FALSE,       # No standardization as requested
  penalty.factor = penalty_factors
)

# Display best tuning parameters
cat("\nBest tuning parameters:\n")
print(enet_model$bestTune)

# -----------------------------------------------
# 6. Model Evaluation
# -----------------------------------------------
cat("\n==== MODEL EVALUATION ====\n")

# Get predictions on the training data
predictions <- predict(enet_model, newdata = X)
predicted_probs <- predict(enet_model, newdata = X, type = "prob")

# Confusion matrix
conf_matrix <- confusionMatrix(
  predictions, 
  factor(y, levels = c(0, 1), labels = c("no", "yes")),
  positive = "yes"
)

# ROC curve
roc_obj <- roc(y, predicted_probs$yes)
auc_value <- auc(roc_obj)

cat("\nConfusion Matrix:\n")
print(conf_matrix)

cat("\nAUC:", round(auc_value, 3), "\n")


# -----------------------------------------------
# 6A. Bootstrap Validation of Model Performance
# -----------------------------------------------
cat("\n==== BOOTSTRAP VALIDATION OF MODEL PERFORMANCE ====\n")

# Number of bootstrap iterations for model validation
n_boot_valid <- 200
boot_performance <- data.frame(
  AUC = numeric(n_boot_valid),
  Accuracy = numeric(n_boot_valid),
  Sensitivity = numeric(n_boot_valid),
  Specificity = numeric(n_boot_valid),
  PPV = numeric(n_boot_valid),
  NPV = numeric(n_boot_valid)
)

set.seed(456)
cat("Performing bootstrap validation...\n")
for(b in 1:n_boot_valid) {
  if(b %% 20 == 0) cat("  Bootstrap iteration:", b, "\n")
  
  # Create bootstrap sample
  pos_idx <- which(y == 1)
  neg_idx <- which(y == 0)
  
  # Sample with replacement within each stratum
  boot_pos <- sample(pos_idx, length(pos_idx), replace = TRUE)
  boot_neg <- sample(neg_idx, length(neg_idx), replace = TRUE)
  
  # Combine and shuffle
  boot_idx <- sample(c(boot_pos, boot_neg))
  
  X_boot <- X[boot_idx, ]
  y_boot <- y[boot_idx]
  
  # Out-of-bag (OOB) samples
  oob_idx <- setdiff(1:nrow(X), unique(boot_idx))
  if(length(oob_idx) == 0) next  # Skip if no OOB samples
  
  X_oob <- X[oob_idx, , drop = FALSE]
  y_oob <- y[oob_idx]
  
  tryCatch({
    # Train model on bootstrap sample
    boot_model <- train(
      x = X_boot,
      y = factor(y_boot, levels = c(0, 1), labels = c("no", "yes")),
      method = "glmnet",
      trControl = trainControl(method = "none"),  # No CV within bootstrap
      tuneGrid = data.frame(
        alpha = enet_model$bestTune$alpha,
        lambda = enet_model$bestTune$lambda
      ),
      family = "binomial",
      standardize = FALSE,
      penalty.factor = penalty_factors
    )
    
    # Predict on OOB samples
    boot_pred_probs <- predict(boot_model, newdata = X_oob, type = "prob")
    boot_pred_class <- predict(boot_model, newdata = X_oob)
    
    # Calculate performance metrics
    boot_roc <- roc(y_oob, boot_pred_probs$yes, quiet = TRUE)
    boot_conf <- confusionMatrix(
      boot_pred_class,
      factor(y_oob, levels = c(0, 1), labels = c("no", "yes")),
      positive = "yes"
    )
    
    # Store results
    boot_performance$AUC[b] <- auc(boot_roc)
    boot_performance$Accuracy[b] <- boot_conf$overall["Accuracy"]
    boot_performance$Sensitivity[b] <- boot_conf$byClass["Sensitivity"]
    boot_performance$Specificity[b] <- boot_conf$byClass["Specificity"]
    boot_performance$PPV[b] <- boot_conf$byClass["Pos Pred Value"]
    boot_performance$NPV[b] <- boot_conf$byClass["Neg Pred Value"]
    
  }, error = function(e) {
    cat("Warning: Error in bootstrap iteration", b, "\n")
    boot_performance[b, ] <- NA
  })
}

# Remove NA rows
boot_performance <- boot_performance[complete.cases(boot_performance), ]

# Calculate bootstrap statistics
boot_stats <- data.frame(
  Metric = names(boot_performance),
  Mean = apply(boot_performance, 2, mean),
  SD = apply(boot_performance, 2, sd),
  CI_Lower = apply(boot_performance, 2, quantile, probs = 0.025),
  CI_Upper = apply(boot_performance, 2, quantile, probs = 0.975)
)

cat("\nBootstrap Validation Results (", nrow(boot_performance), " successful iterations):\n")
print(boot_stats)

# Save bootstrap validation results
write.csv(boot_stats, "CMML_Bootstrap_Model_Performance.csv", row.names = FALSE)

# -----------------------------------------------
# 6B. Optimism-Corrected Internal Validation
# -----------------------------------------------
cat("\n==== OPTIMISM-CORRECTED INTERNAL VALIDATION ====\n")

# Calculate apparent performance (on full dataset)
apparent_pred_probs <- predict(enet_model, newdata = X, type = "prob")
apparent_roc <- roc(y, apparent_pred_probs$yes, quiet = TRUE)
apparent_auc <- auc(apparent_roc)

# Optimism-corrected bootstrap
n_optimism_boot <- 200
optimism_values <- numeric(n_optimism_boot)

set.seed(789)
cat("Calculating optimism correction...\n")
for(b in 1:n_optimism_boot) {
  if(b %% 20 == 0) cat("  Bootstrap iteration:", b, "\n")
  
  # Create bootstrap sample
  pos_idx <- which(y == 1)
  neg_idx <- which(y == 0)
  
  # Sample with replacement within each stratum
  boot_pos <- sample(pos_idx, length(pos_idx), replace = TRUE)
  boot_neg <- sample(neg_idx, length(neg_idx), replace = TRUE)
  
  # Combine and shuffle
  boot_idx <- sample(c(boot_pos, boot_neg))
  
  X_boot <- X[boot_idx, ]
  y_boot <- y[boot_idx]
  
  tryCatch({
    # Train on bootstrap sample
    boot_model <- train(
      x = X_boot,
      y = factor(y_boot, levels = c(0, 1), labels = c("no", "yes")),
      method = "glmnet",
      trControl = trainControl(method = "none"),
      tuneGrid = data.frame(
        alpha = enet_model$bestTune$alpha,
        lambda = enet_model$bestTune$lambda
      ),
      family = "binomial",
      standardize = FALSE,
      penalty.factor = penalty_factors
    )
    
    # Performance on bootstrap sample (apparent in bootstrap)
    boot_pred_on_boot <- predict(boot_model, newdata = X_boot, type = "prob")
    boot_auc_on_boot <- auc(roc(y_boot, boot_pred_on_boot$yes, quiet = TRUE))
    
    # Performance on original sample
    boot_pred_on_orig <- predict(boot_model, newdata = X, type = "prob")
    boot_auc_on_orig <- auc(roc(y, boot_pred_on_orig$yes, quiet = TRUE))
    
    # Optimism = performance on bootstrap - performance on original
    optimism_values[b] <- boot_auc_on_boot - boot_auc_on_orig
    
  }, error = function(e) {
    optimism_values[b] <- NA
  })
}

# Calculate optimism-corrected performance
mean_optimism <- mean(optimism_values, na.rm = TRUE)
optimism_corrected_auc <- apparent_auc - mean_optimism

cat("\nOptimism-Corrected Performance:\n")
cat("Apparent AUC:", round(apparent_auc, 3), "\n")
cat("Mean optimism:", round(mean_optimism, 3), "\n")
cat("Optimism-corrected AUC:", round(optimism_corrected_auc, 3), "\n")
cat("95% CI for optimism:", round(quantile(optimism_values, c(0.025, 0.975), na.rm = TRUE), 3), "\n")

# -----------------------------------------------
# 6C. Repeated Cross-Validation for Stability
# -----------------------------------------------
cat("\n==== REPEATED CROSS-VALIDATION FOR MODEL STABILITY ====\n")

# Perform repeated CV
set.seed(321)
repeated_cv_control <- trainControl(
  method = "repeatedcv",
  number = 10,              # 10-fold CV
  repeats = 10,            # 10 repeats
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

cat("Performing 10x10 repeated cross-validation...\n")
repeated_cv_model <- train(
  x = X,
  y = factor(y, levels = c(0, 1), labels = c("no", "yes")),
  method = "glmnet",
  trControl = repeated_cv_control,
  tuneGrid = data.frame(
    alpha = enet_model$bestTune$alpha,
    lambda = enet_model$bestTune$lambda
  ),
  metric = "ROC",
  family = "binomial",
  standardize = FALSE,
  penalty.factor = penalty_factors
)

# Extract CV results
cv_results <- repeated_cv_model$resample
cv_summary <- data.frame(
  Metric = c("ROC/AUC", "Sensitivity", "Specificity"),
  Mean = c(mean(cv_results$ROC), mean(cv_results$Sens), mean(cv_results$Spec)),
  SD = c(sd(cv_results$ROC), sd(cv_results$Sens), sd(cv_results$Spec)),
  Min = c(min(cv_results$ROC), min(cv_results$Sens), min(cv_results$Spec)),
  Max = c(max(cv_results$ROC), max(cv_results$Sens), max(cv_results$Spec))
)

cat("\nRepeated Cross-Validation Results:\n")
print(cv_summary)

# Create visualization of validation results
pdf("CMML_Model_Validation_Results.pdf", width = 12, height = 8)
par(mfrow = c(2, 2))

# 1. Bootstrap AUC distribution
hist(boot_performance$AUC, breaks = 20, main = "Bootstrap AUC Distribution",
     xlab = "AUC", col = "lightblue", border = "darkblue")
abline(v = mean(boot_performance$AUC), col = "red", lwd = 2, lty = 2)
abline(v = quantile(boot_performance$AUC, c(0.025, 0.975)), col = "red", lwd = 1, lty = 3)

# 2. Bootstrap performance metrics
boxplot(boot_performance[, c("Sensitivity", "Specificity", "PPV", "NPV")],
        main = "Bootstrap Performance Metrics", ylab = "Value",
        col = c("lightgreen", "lightcoral", "lightyellow", "lightpink"))

# 3. Optimism distribution
hist(optimism_values[!is.na(optimism_values)], breaks = 20,
     main = "Optimism Distribution", xlab = "Optimism",
     col = "lightgray", border = "darkgray")
abline(v = mean_optimism, col = "red", lwd = 2, lty = 2)

# 4. Repeated CV ROC distribution
hist(cv_results$ROC, breaks = 20, main = "Repeated CV AUC Distribution",
     xlab = "AUC", col = "lightcyan", border = "darkblue")
abline(v = mean(cv_results$ROC), col = "red", lwd = 2, lty = 2)

dev.off()

# Save all validation results
validation_summary <- list(
  Bootstrap_Performance = boot_stats,
  Optimism_Corrected = data.frame(
    Apparent_AUC = apparent_auc,
    Mean_Optimism = mean_optimism,
    Corrected_AUC = optimism_corrected_auc
  ),
  Repeated_CV = cv_summary
)

saveRDS(validation_summary, "CMML_Model_Validation_Summary.rds")
cat("\nAll validation results saved to CMML_Model_Validation_Summary.rds\n")
cat("\nValidation plots saved to CMML_Model_Validation_Results.pdf\n")

lapply(validation_summary, function(x) write.table( data.frame(x), 'test.csv'  , append= T, sep=',' ))



# -----------------------------------------------
# 7. Extract Model Coefficients
# -----------------------------------------------
cat("\n==== MODEL COEFFICIENTS ====\n")

# Extract coefficients from the final model
coefs <- coef(enet_model$finalModel, s = enet_model$bestTune$lambda)
coefs_df <- as.data.frame(as.matrix(coefs))
colnames(coefs_df) <- "Coefficient"
coefs_df$Variable <- rownames(coefs)

# Sort by absolute coefficient value
coefs_df <- coefs_df[order(abs(coefs_df$Coefficient), decreasing = TRUE), ]

# Display non-zero coefficients
non_zero_coefs <- coefs_df[coefs_df$Coefficient != 0, ]
cat("\nNon-zero coefficients:\n")
print(non_zero_coefs)

# Calculate odds ratios
non_zero_coefs$OddsRatio <- exp(non_zero_coefs$Coefficient)

# Save to CSV

# -----------------------------------------------
# 9. Create a Clinical Prediction Model
# -----------------------------------------------
cat("\n==== DEVELOPING CLINICAL PREDICTION MODEL ====\n")

# Create a point-based scoring system from non-zero coefficient predictors
model_vars <- non_zero_coefs[non_zero_coefs$Variable != "(Intercept)", ]

if (nrow(model_vars) > 0 && any(model_vars$Coefficient != 0)) {
  # Ensure there's at least one variable with a non-zero coefficient
  # Calculate points by scaling coefficients
  # The variable with the largest absolute coefficient effect gets approx. +5 or -5 points
  max_abs_coef <- max(abs(model_vars$Coefficient[model_vars$Coefficient != 0])) # Max of non-zero absolute coefficients
  
  if (is.finite(max_abs_coef) && max_abs_coef > 0) {
    model_vars$Points <- round(5 * model_vars$Coefficient / max_abs_coef)
    # Handle cases where some coefficients might be zero if non_zero_coefs had very small values rounded by R's printing # nolint # nolint: line_length_linter.
    model_vars$Points[model_vars$Coefficient == 0] <- 0 
  } else {
    # Fallback if max_abs_coef is not suitable for scaling (e.g., all remaining are zero, which shouldn't happen here)
    warning("Could not scale points appropriately. Max absolute coefficient was zero or non-finite.")
    model_vars$Points <- 0 
  }

  # Print the clinical prediction model (points table)
  cat("\nClinical Prediction Model for CMML Progression (Risk Factors & Points):\n")
  cat("---------------------------------------------------------------------\n")
  cat(sprintf("%-40s | %s\n", "Risk Factor", "Points"))
  cat("---------------------------------------------------------------------\n")
  for (i in 1:nrow(model_vars)) {
    cat(sprintf("%-40s | %+d\n", model_vars$Variable[i], model_vars$Points[i]))
  }
  cat("---------------------------------------------------------------------\n")

  # Save the clinical prediction model variables and points
  write.csv(model_vars[, c("Variable", "Coefficient", "OddsRatio", "Points")],
            "CMML_Clinical_Prediction_Model_Variables.csv", row.names = FALSE)
  cat("\nPoint system (variables and points) saved to CMML_Clinical_Prediction_Model_Variables.csv\n")

  ###############################################
  ### START OF REVISED CODE BLOCK             ###
  ###############################################
  # (The new code provided below will replace this comment line and the one after it)
  ### TODO: add prediction for clinicians     ###
  ###############################################

  # (The rest of your script, e.g., section 10. Bootstrap, follows after this block)

} else {
  cat("\nNo predictors (excluding intercept) selected by the Elastic Net model.\n")
  cat("Cannot create a point-based scoring system or identify a total score cut-off.\n")
}


cat("\n==== IDENTIFYING AND EVALUATING TOTAL SCORE CUT-OFF FOR CLINICAL PREDICTION ====\n")

  # 1. Calculate total scores for each patient
  # Ensure 'X' (predictors_imputed matrix) and 'y' (outcome vector) are the ones used for model training.
  # 'model_vars' contains the variables and their points.

  predictor_names_in_point_model <- model_vars$Variable
  
  # Select the corresponding columns from X. Ensure X has these columns.
  # It's crucial that X here is the same matrix used for training the enet_model
  # X <- as.matrix(predictors_imputed) # This should have been defined before model training
  
  # Check if all variables in the point model are present in X
  if (!all(predictor_names_in_point_model %in% colnames(X))) {
    stop("Mismatch: Not all variables from the point model are found in the predictor matrix X. Check variable names.")
  }
  X_for_scoring <- X[, predictor_names_in_point_model, drop = FALSE]
  
  # Get the points for these predictors
  points_for_scoring <- model_vars$Points
  
  # Calculate total score for each patient: (patient_1_var1 * points_var1) + (patient_1_var2 * points_var2) + ...
  total_scores <- as.vector(X_for_scoring %*% points_for_scoring)
  
  cat("Summary of calculated total scores for patients:\n")
  print(summary(total_scores))
  # Visualize distribution of total scores
  pdf("CMML_Total_Scores_Distribution.pdf", width = 8, height = 6)
  hist(total_scores, 
       main = "Distribution of Patient Total Risk Scores", 
       xlab = "Total Score", 
       border = "blue", col = "lightblue", breaks = 20)
  dev.off()
  cat("\nDistribution plot of total scores saved to CMML_Total_Scores_Distribution.pdf\n")

  # 2. Find Optimal Cut-off for the Total Score
  # Ensure 'y' (outcome) is available (0 for no progression, 1 for progression)
  scoring_cutpoint_data <- data.frame(
    TotalScore = total_scores,
    Progression = factor(y, levels = c(0, 1), labels = c("no", "yes")) # 'y' is the original outcome variable
  )
  
  # Check for variability in scores
  if (length(unique(scoring_cutpoint_data$TotalScore)) > 1) {
    cat("\nDetermining optimal cut-point for total score using OptimalCutpoints package...\n")
    # Find optimal cutpoint using multiple methods for comparison, focusing on Youden
    score_cut_result <- OptimalCutpoints::optimal.cutpoints(
      X = "TotalScore",
      status = "Progression",
      tag.healthy = "no", # "no" corresponds to outcome 0
      methods = c("Youden", "SpEqualSe", "MaxSp", "MaxSe", "PROC01"), # PROC01 minimizes distance to (0,1) on ROC
      data = scoring_cutpoint_data,
      ci.fit = TRUE, 
      conf.level = 0.95,
      control = control.cutpoints(CFP = 0) # Control for Cost Function: assume equal misclassification costs here, or adjust
    )
    
    cat("\nSummary of Optimal Total Score Cut-offs (from OptimalCutpoints):\n")
    print(summary(score_cut_result))
    
    # Select a primary cut-off, e.g., Youden's J index
    if (!is.null(score_cut_result$Youden) && 
        !is.null(score_cut_result$Youden$Global) &&
        !is.null(score_cut_result$Youden$Global$optimal.cutoff) &&
        length(score_cut_result$Youden$Global$optimal.cutoff$cutoff) > 0) {
          
      optimal_total_score_cutoff <- score_cut_result$Youden$Global$optimal.cutoff$cutoff[1] # Take the first if multiple
      cat("\nSelected Optimal Total Score Cut-off (Youden Index):", round(optimal_total_score_cutoff, 2), "\n")
      
      # 3. Evaluate the performance of this chosen cut-off
      predicted_progression_by_score <- ifelse(total_scores >= optimal_total_score_cutoff, "yes", "no")
      predicted_progression_by_score <- factor(predicted_progression_by_score, levels = c("no", "yes"))
      actual_progression <- factor(y, levels = c(0, 1), labels = c("no", "yes"))
      
      score_conf_matrix <- confusionMatrix(
        data = predicted_progression_by_score,
        reference = actual_progression,
        positive = "yes"
      )
      cat("\nConfusion Matrix for Total Score-Based Prediction (Cut-off =", round(optimal_total_score_cutoff, 2), "):\n")
      print(score_conf_matrix)
      
      # ROC curve for the total score model itself
      score_roc_obj <- pROC::roc(response = actual_progression, predictor = total_scores, quiet = TRUE, levels=c("no", "yes"), direction="<")
      score_auc_value <- pROC::auc(score_roc_obj)
      cat("\nAUC for the Total Score model:", round(score_auc_value, 3), "\n")
      
      # Plot ROC curve for the total score
      pdf("CMML_TotalScore_ROC_Curve.pdf", width = 8, height = 8)
      plot(score_roc_obj, main = "ROC Curve for CMML Progression (Total Score Model)",
           col = "darkgreen", lwd = 2, print.thres = "best", print.thres.best.method="youden")
      #abline(a = 0, b = 1, lty = 2, col = "gray") # Line of no discrimination
      legend("bottomright", legend = paste("AUC =", round(score_auc_value, 3)), bty = "n")
      # Highlighting the specific Youden cut-off chosen from OptimalCutpoints
      youden_sens <- score_cut_result$Youden$Global$optimal.cutoff$Se[[1]]
      youden_spec <- score_cut_result$Youden$Global$optimal.cutoff$Sp[[1]]
      points(x = 1 - youden_spec, y = youden_sens, pch = 19, cex = 2, col = "red")
      text(x = 1 - youden_spec, y = youden_sens, 
           labels = paste0("Cutoff: ", round(optimal_total_score_cutoff, 2),
                          "\nSens: ", round(youden_sens, 3), 
                          "\nSpec: ", round(youden_spec, 3)),
           pos = 4, cex = 0.9, col = "red")
      dev.off()
      cat("\nROC curve for total score model saved to CMML_TotalScore_ROC_Curve.pdf\n")
      
      # 4. Clinical Application Guidance
      cat("\n---------------------------------------------------------------------\n")
      cat("           CLINICAL APPLICATION GUIDANCE (Point System)\n")
      cat("---------------------------------------------------------------------\n")
      cat("1. For each patient, calculate the Total Risk Score using the point system table.\n")
      cat(paste0("2. If Total Score >= ", round(optimal_total_score_cutoff, 2), ": Predict HIGH RISK for CMML progression.\n"))
      cat(paste0("3. If Total Score < ", round(optimal_total_score_cutoff, 2), ": Predict LOW RISK for CMML progression.\n"))
      cat("\nPerformance of this rule (based on Youden Index cut-off):\n")
      cat(sprintf("  - Sensitivity: %.3f\n", score_conf_matrix$byClass["Sensitivity"]))
      cat(sprintf("  - Specificity: %.3f\n", score_conf_matrix$byClass["Specificity"]))
      cat(sprintf("  - Positive Predictive Value (PPV): %.3f\n", score_conf_matrix$byClass["Pos Pred Value"]))
      cat(sprintf("  - Negative Predictive Value (NPV): %.3f\n", score_conf_matrix$byClass["Neg Pred Value"]))
      cat(sprintf("  - Accuracy: %.3f\n", score_conf_matrix$overall["Accuracy"]))
      cat("---------------------------------------------------------------------\n")
      
      # Optionally, save the cutoff value and its performance
      cutoff_summary <- data.frame(
          Cutoff_Value = optimal_total_score_cutoff,
          Method = "Youden Index",
          AUC_TotalScore = score_auc_value,
          Sensitivity = score_conf_matrix$byClass["Sensitivity"],
          Specificity = score_conf_matrix$byClass["Specificity"],
          PPV = score_conf_matrix$byClass["Pos Pred Value"],
          NPV = score_conf_matrix$byClass["Neg Pred Value"],
          Accuracy = score_conf_matrix$overall["Accuracy"]
      )
      write.csv(cutoff_summary, "CMML_TotalScore_Cutoff_Performance.csv", row.names = FALSE)
      cat("\nTotal score cut-off and performance summary saved to CMML_TotalScore_Cutoff_Performance.csv\n")
      
    } else {
      cat("\nCould not determine an optimal cut-off for the total score using Youden Index (e.g., method failed or no valid cut-off found).\n")
      cat("Please examine the distribution of scores and the output from OptimalCutpoints.\n")
      cat("You might need to select a cut-off based on other criteria or methods if Youden is not appropriate.\n")
    }
  } else {
    cat("\nNot enough variability in total scores (all patients have the same or very few unique scores).\n")
    cat("Cannot determine a meaningful cut-off. Distribution of total scores:\n")
    print(table(total_scores))
  }





# ===============================================
# PUBLICATION-QUALITY FIGURE FUNCTIONS AND SETUP
# ===============================================

# Load additional libraries for publication figures
library(ggplot2)
library(cowplot)
library(scales)
library(viridis)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)

# Define publication theme
theme_publication <- function(base_size = 12, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # Panel
      panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
      panel.grid.major = element_line(colour = "grey90", size = 0.25),
      panel.grid.minor = element_blank(),
      
      # Axes
      axis.text = element_text(size = rel(0.9), colour = "black"),
      axis.title = element_text(size = rel(1.0), face = "bold"),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.line = element_line(colour = "black", size = 0.5),
      
      # Legend
      legend.key = element_rect(fill = "white", colour = NA),
      legend.title = element_text(size = rel(0.9), face = "bold"),
      legend.text = element_text(size = rel(0.8)),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key.size = unit(1.2, "lines"),
      
      # Strips for facets
      strip.background = element_rect(fill = "grey95", colour = "black", size = 0.5),
      strip.text = element_text(size = rel(0.9), face = "bold"),
      
      # Overall
      plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0.5,
                               margin = margin(b = 10)),
      plot.subtitle = element_text(size = rel(1.0), hjust = 0.5,
                                  margin = margin(b = 10)),
      plot.margin = unit(c(1, 1, 1, 1), "lines")
    )
}

# Color palette for publication
pub_colors <- c(
  "positive" = "#D73027",  # Red for positive coefficients
  "negative" = "#4575B4",  # Blue for negative coefficients
  "neutral" = "#969696",   # Gray
  "highlight" = "#FF7F00", # Orange for emphasis
  "significant" = "#1B9E77", # Teal for significant
  "nonsignificant" = "#CCCCCC" # Light gray
)


# ===============================================
# REPLACE SECTION 7.3: ROC CURVE
# ===============================================

# Calculate ROC curve data
roc_data <- data.frame(
  sensitivity = roc_obj$sensitivities,
  specificity = roc_obj$specificities,
  fpr = 1 - roc_obj$specificities
)

# Find optimal point (Youden)
optimal_idx <- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
optimal_point <- data.frame(
  fpr = 1 - roc_obj$specificities[optimal_idx],
  sensitivity = roc_obj$sensitivities[optimal_idx],
  threshold = roc_obj$thresholds[optimal_idx]
)

# Figure 3: ROC Curve
fig3_roc <- ggplot(roc_data, aes(x = fpr, y = sensitivity)) +
  geom_line(size = 1.2, color = pub_colors["positive"]) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "gray50", size = 0.8) +
  geom_point(data = optimal_point, aes(x = fpr, y = sensitivity),
             size = 4, shape = 21, fill = pub_colors["highlight"],
             color = "black", stroke = 1.5) +
  annotate("text", x = 0.6, y = 0.3,
           label = paste0("AUC = ", sprintf("%.3f", auc_value)),
           size = 5, fontface = "bold") +
  annotate("text", x = optimal_point$fpr + 0.05, y = optimal_point$sensitivity - 0.05,
           label = paste0("Optimal cutoff\n",
                         "Sensitivity: ", sprintf("%.2f", optimal_point$sensitivity), "\n",
                         "Specificity: ", sprintf("%.2f", 1 - optimal_point$fpr)),
           size = 3.5, hjust = 0) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
  labs(title = "ROC Curve for CMML Progression Prediction Model",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_publication() +
  coord_equal()

ggsave("Figure3_CMML_ROC.pdf", fig3_roc,
       width = 8, height = 8, dpi = 300, device = cairo_pdf)
ggsave("Figure3_CMML_ROC.tiff", fig3_roc,
       width = 8, height = 8, dpi = 300, compression = "lzw")

# ===============================================
# REPLACE VALIDATION RESULTS PLOTS
# ===============================================

if(exists("boot_performance")) {
  # Create a comprehensive validation figure with multiple panels
  
  # Panel A: Bootstrap AUC Distribution
  panel_a <- ggplot(boot_performance, aes(x = AUC)) +
    geom_histogram(aes(y = ..density..), bins = 25,
                   fill = pub_colors["positive"], alpha = 0.7,
                   color = "black", size = 0.3) +
    geom_density(size = 1, color = "black") +
    geom_vline(xintercept = mean(boot_performance$AUC),
               linetype = "dashed", size = 1, color = "black") +
    geom_vline(xintercept = quantile(boot_performance$AUC, c(0.025, 0.975)),
               linetype = "dotted", size = 0.8, color = "black") +
    labs(title = "A. Bootstrap AUC Distribution",
         x = "AUC",
         y = "Density") +
    theme_publication(base_size = 10) +
    annotate("text", x = mean(boot_performance$AUC),
             y = max(density(boot_performance$AUC)$y) * 0.9,
             label = sprintf("Mean = %.3f", mean(boot_performance$AUC)),
             hjust = -0.1, size = 3.5)
  
  # Panel B: Performance Metrics Comparison
  perf_metrics_long <- tidyr::pivot_longer(
    boot_performance[, c("Sensitivity", "Specificity", "PPV", "NPV")],
    cols = everything(),
    names_to = "Metric",
    values_to = "Value"
  )
  
  panel_b <- ggplot(perf_metrics_long, aes(x = Metric, y = Value, fill = Metric)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "B. Bootstrap Performance Metrics",
         x = NULL,
         y = "Value") +
    theme_publication(base_size = 10)
  
  # Panel C: Calibration Plot (if predictions available)
  if(exists("predicted_probs")) {
    # Create calibration data
    cal_data <- data.frame(
      pred_prob = predicted_probs$yes,
      actual = y
    )
    cal_data$bin <- cut(cal_data$pred_prob,
                       breaks = seq(0, 1, 0.1),
                       include.lowest = TRUE)
    
    cal_summary <- cal_data %>%
      group_by(bin) %>%
      summarise(
        mean_pred = mean(pred_prob),
        mean_actual = mean(actual),
        n = n(),
        se = sqrt(mean_actual * (1 - mean_actual) / n)
      )
    
    panel_c <- ggplot(cal_summary, aes(x = mean_pred, y = mean_actual)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                  color = "gray50", size = 1) +
      geom_errorbar(aes(ymin = mean_actual - 1.96*se,
                       ymax = mean_actual + 1.96*se),
                    width = 0.02, color = "black") +
      geom_point(aes(size = n), color = pub_colors["positive"], alpha = 0.8) +
      scale_size_continuous(range = c(3, 8), name = "n") +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(title = "C. Calibration Plot",
           x = "Predicted Probability",
           y = "Observed Proportion") +
      theme_publication(base_size = 10) +
      coord_equal()
  }
  
  # Panel D: Optimism Distribution (if available)
  if(exists("optimism_values")) {
    panel_d <- ggplot(data.frame(optimism = optimism_values[!is.na(optimism_values)]),
                      aes(x = optimism)) +
      geom_histogram(aes(y = ..density..), bins = 25,
                     fill = pub_colors["neutral"], alpha = 0.7,
                     color = "black", size = 0.3) +
      geom_density(size = 1, color = "black") +
      geom_vline(xintercept = mean_optimism,
                 linetype = "dashed", size = 1, color = "black") +
      labs(title = "D. Optimism Distribution",
           x = "Optimism",
           y = "Density") +
      theme_publication(base_size = 10) +
      annotate("text", x = mean_optimism,
               y = max(density(optimism_values[!is.na(optimism_values)])$y) * 0.9,
               label = sprintf("Mean = %.3f", mean_optimism),
               hjust = -0.1, size = 3.5)
  }
  
  # Combine panels
  if(exists("panel_c") && exists("panel_d")) {
    fig6_validation <- cowplot::plot_grid(panel_a, panel_b, panel_c, panel_d,
                                         ncol = 2, nrow = 2,
                                         align = "hv")
  } else {
    fig6_validation <- cowplot::plot_grid(panel_a, panel_b,
                                         ncol = 2, nrow = 1,
                                         align = "hv")
  }
  
  # Add overall title
  fig6_validation_final <- cowplot::ggdraw() +
    cowplot::draw_plot(fig6_validation, 0, 0, 1, 0.95) +
    cowplot::draw_label("Model Validation Results",
                       x = 0.5, y = 0.98,
                       size = 14, fontface = "bold",
                       hjust = 0.5)
  
  ggsave("Figure6_CMML_Validation.pdf", fig6_validation_final,
         width = 12, height = 10, dpi = 300, device = cairo_pdf)
  ggsave("Figure6_CMML_Validation.tiff", fig6_validation_final,
         width = 12, height = 10, dpi = 300, compression = "lzw")
}

# ===============================================
# CREATE SUPPLEMENTARY FIGURES
# ===============================================

# Supplementary Figure 1: Variable Importance (if available)
if(exists("enet_model")) {
  # Extract variable importance
  var_imp <- varImp(enet_model, scale = FALSE)
  var_imp_data <- data.frame(
    Variable = rownames(var_imp$importance),
    Importance = var_imp$importance[,1]
  )
  var_imp_data <- var_imp_data[order(var_imp_data$Importance, decreasing = TRUE), ]
  var_imp_data <- head(var_imp_data, 30)
  
  supp_fig1 <- ggplot(var_imp_data,
                      aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_segment(aes(xend = Variable, yend = 0), size = 1, color = "gray50") +
    geom_point(size = 3, color = pub_colors["positive"]) +
    coord_flip() +
    labs(title = "Variable Importance in CMML Progression Model",
         subtitle = "Top 30 variables by absolute coefficient value",
         x = NULL,
         y = "Importance (Absolute Coefficient)") +
    theme_publication()
  
  ggsave("Supplementary_Figure1_VarImportance.pdf", supp_fig1,
         width = 8, height = 10, dpi = 300, device = cairo_pdf)
}

# Supplementary Figure 2: Cross-validation performance across alpha values
if(exists("enet_model") && !is.null(enet_model$results)) {
  cv_results_plot <- enet_model$results %>%
    group_by(alpha) %>%
    summarise(
      mean_AUC = mean(ROC),
      sd_AUC = sd(ROC),
      .groups = "drop"
    )
  
  supp_fig2 <- ggplot(cv_results_plot, aes(x = alpha, y = mean_AUC)) +
    geom_ribbon(aes(ymin = mean_AUC - sd_AUC, ymax = mean_AUC + sd_AUC),
                alpha = 0.2, fill = pub_colors["positive"]) +
    geom_line(size = 1.2, color = pub_colors["positive"]) +
    geom_point(size = 3, color = pub_colors["positive"]) +
    geom_vline(xintercept = enet_model$bestTune$alpha,
               linetype = "dashed", color = "black", size = 1) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    labs(title = "Cross-Validation Performance Across Alpha Values",
         subtitle = "Alpha = 0 (Ridge) to Alpha = 1 (Lasso)",
         x = "Alpha (Elastic Net Mixing Parameter)",
         y = "Mean AUC ± SD") +
    theme_publication() +
    annotate("text", x = enet_model$bestTune$alpha,
             y = min(cv_results_plot$mean_AUC) * 0.98,
             label = paste0("Optimal α = ", enet_model$bestTune$alpha),
             hjust = -0.1, fontface = "bold")
  
  ggsave("Supplementary_Figure2_AlphaSelection.pdf", supp_fig2,
         width = 8, height = 6, dpi = 300, device = cairo_pdf)
}

cat("\nAll publication-quality figures have been generated in PDF and TIFF formats.\n")
cat("Figures are saved at 300 DPI for journal submission requirements.\n")





# -----------------------------------------------
# 8. Visualizing Results
# -----------------------------------------------
cat("\n==== CREATING VISUALIZATIONS ====\n")

# Create coefficient plot
pdf("CMML_ElasticNet_Coefficients.pdf", width = 10, height = 8)
ggplot(head(non_zero_coefs[non_zero_coefs$Variable != "(Intercept)", ], 20), 
       aes(x = reorder(Variable, abs(Coefficient)), y = Coefficient)) +
  geom_bar(stat = "identity", aes(fill = Coefficient > 0)) +
  scale_fill_manual(values = c("firebrick", "dodgerblue"), 
                    labels = c("Negative", "Positive"),
                    name = "Direction") +
  coord_flip() +
  labs(title = "Elastic Net - Top Variables by Coefficient Magnitude",
       x = "Variable",
       y = "Coefficient Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Create odds ratio plot
ggplot(head(non_zero_coefs[non_zero_coefs$Variable != "(Intercept)", ], 20), 
       aes(x = reorder(Variable, OddsRatio), y = OddsRatio)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_log10() +
  coord_flip() +
  labs(title = "Odds Ratios for CMML Progression Predictors",
       x = "Variable",
       y = "Odds Ratio (log scale)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
dev.off()

# ROC curve
pdf("CMML_ElasticNet_ROC.pdf", width = 8, height = 8)
plot(roc_obj, main = "ROC Curve for CMML Progression Model",
     col = "blue", lwd = 2)
text(0.6, 0.2, paste("AUC =", round(auc_value, 3)), cex = 1.2)
abline(a = 1, b = -1, lty = 2, col = "gray")
dev.off()

# -----------------------------------------------
# 9. Create a Clinical Prediction Model
# -----------------------------------------------
cat("\n==== DEVELOPING CLINICAL PREDICTION MODEL ====\n")
## TODO: help me Identify a total score cut-off that allows to accurately predict CMML progression so that we can include it as part of the more clinically applicable model using points.

# Create a point-based scoring system
# Assign points based on coefficient magnitude
model_vars <- non_zero_coefs[non_zero_coefs$Variable != "(Intercept)", ]
max_coef <- max(abs(model_vars$Coefficient))
model_vars$Points <- round(5 * model_vars$Coefficient / max_coef)

# Print the clinical prediction model
cat("\nClinical Prediction Model for CMML Progression:\n")
cat("---------------------------------------------\n")
cat("Risk Factor                      | Points\n")
cat("---------------------------------------------\n")
for (i in 1:nrow(model_vars)) {
  cat(sprintf("%-30s | %+d\n", model_vars$Variable[i], model_vars$Points[i]))
}
cat("---------------------------------------------\n")

# Save the clinical prediction model
write.csv(model_vars[, c("Variable", "Coefficient", "OddsRatio", "Points")],
          "CMML_Clinical_Prediction_Model.csv", row.names = FALSE)

###############################################
### TODO: add prediction for clinicians     ###
###############################################

# -----------------------------------------------
# 10. Bootstrap for Confidence Intervals
# -----------------------------------------------
cat("\n==== CALCULATING BOOTSTRAP CONFIDENCE INTERVALS ====\n")

# Number of bootstrap iterations
n_bootstrap <- 500
bootstrap_results <- matrix(NA, nrow = n_bootstrap, ncol = ncol(X) + 1)  # +1 for intercept
colnames(bootstrap_results) <- c("(Intercept)", colnames(X))

set.seed(42)
for(i in 1:n_bootstrap) {
  if(i %% 50 == 0) cat("  Bootstrap iteration:", i, "\n")
  
  # Sample with replacement
  boot_idx <- sample(nrow(X), replace = TRUE)
  X_boot <- X[boot_idx, ]
  y_boot <- y[boot_idx]
  
  # Fit Elastic Net with fixed parameters from the original model
  tryCatch({
    boot_model <- glmnet(
      x = X_boot,
      y = y_boot,
      alpha = enet_model$bestTune$alpha,
      lambda = enet_model$bestTune$lambda,
      family = "binomial"
    )
    
    # Extract coefficients
    coef_boot <- as.vector(coef(boot_model))
    bootstrap_results[i, ] <- coef_boot
  }, error = function(e) {
    cat("Warning: Error in bootstrap iteration", i, ":", conditionMessage(e), "\n")
  })
}

# Calculate 95% confidence intervals
ci_lower <- apply(bootstrap_results, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
ci_upper <- apply(bootstrap_results, 2, function(x) quantile(x, 0.975, na.rm = TRUE))

# Calculate odds ratios and CI's
or_lower <- exp(ci_lower)
or_upper <- exp(ci_upper)

# Create bootstrap results data frame
bootstrap_inference <- data.frame(
  Variable = names(ci_lower),
  Coefficient = as.vector(coef(enet_model$finalModel, 
                              s = enet_model$bestTune$lambda)),
  CI_Lower = ci_lower,
  CI_Upper = ci_upper,
  OddsRatio = exp(as.vector(coef(enet_model$finalModel, 
                                s = enet_model$bestTune$lambda))),
  OR_CI_Lower = or_lower,
  OR_CI_Upper = or_upper
)

# Determine significance based on confidence intervals
bootstrap_inference$Significant <- (bootstrap_inference$CI_Lower > 0 & 
                                   bootstrap_inference$CI_Upper > 0) | 
                                  (bootstrap_inference$CI_Lower < 0 & 
                                   bootstrap_inference$CI_Upper < 0)

# Filter to significant variables
significant_vars <- bootstrap_inference[bootstrap_inference$Significant & 
                                       bootstrap_inference$Variable != "(Intercept)", ]
significant_vars <- significant_vars[order(abs(significant_vars$Coefficient), 
                                          decreasing = TRUE), ]

# Print results for significant variables
cat("\nSignificant predictors with confidence intervals:\n")
print(significant_vars[, c("Variable", "Coefficient", "OddsRatio", 
                         "OR_CI_Lower", "OR_CI_Upper")])

# Save results
write.csv(significant_vars, "CMML_Progression_Significant_Predictors.csv", row.names = FALSE)

# Create plot with confidence intervals
pdf("CMML_Bootstrap_Inference.pdf", width = 10, height = 8)
ggplot(head(significant_vars, 20), aes(x = reorder(Variable, abs(Coefficient)), y = Coefficient)) +
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

cat("\nAnalysis completed.\n")
