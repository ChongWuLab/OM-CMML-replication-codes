
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

data2 <- read.csv("OM-CMML for time-to-event penalized regression.csv", stringsAsFactors = FALSE)

data2 <- data2[, c("PatientID","Time.to.event.CMML")]


library(dplyr)

# Inner join: keep only patients present in BOTH data1 and data2
combined <- data %>%
  inner_join(data2, by = c("Patient_ID" = "PatientID"))

survtime = combined$Time.to.event.CMML


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
survtime = survtime[!indx]
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

library(survival)
library(glmnet)

# Prepare data for modeling
X <- as.matrix(predictors_imputed)

y_surv <- Surv(time = survtime, event = outcome)

# -------------------------------------------------
fit_cox <- glmnet(
  x      = X,
  y      = y_surv,
  family = "cox",    # tells glmnet to fit Cox model
  alpha  = 0.5       # elastic net; 1 = lasso, 0 = ridge
)

# `fit_cox` contains the full regularization path (the “cox path”)
print(fit_cox)
plot(fit_cox)

# -------------------------------------------------
# 2) Cross-validated Cox elastic net
# -------------------------------------------------
cv_cox <- cv.glmnet(
  x      = X,
  y      = y_surv,
  family = "cox",
  alpha  = 0.5,
  nfolds = 10
)

# Best lambdas
cv_cox$lambda.min      # lambda with minimum CV error
cv_cox$lambda.1se      # more regularized choice

# Coefficients at selected lambda
coef_min  <- coef(cv_cox, s = "lambda.min")
coef_1se  <- coef(cv_cox, s = "lambda.1se")



# Coefficients at lambda.min
beta_min <- coef(cv_cox, s = "lambda.min")
coef_vec <- as.numeric(beta_min)
var_names <- rownames(beta_min)

nonzero_idx  <- which(coef_vec != 0)
nonzero_vars <- var_names[nonzero_idx]

coefs_df <- data.frame(
  Variable    = var_names,
  Coefficient = coef_vec
)
coefs_df <- coefs_df[coefs_df$Coefficient != 0 & coefs_df$Variable != "(Intercept)", ]
coefs_df$HazardRatio <- exp(coefs_df$Coefficient)
coefs_df <- coefs_df[order(abs(coefs_df$Coefficient), decreasing = TRUE), ]

cat("\nNon-zero Cox coefficients (lambda.min):\n")
print(head(coefs_df, 30))

write.csv(
  coefs_df,
  file = "CMML_Cox_Clinical_Prediction_Model_Variables.csv",
  row.names = FALSE
)


library(Hmisc)

# Linear predictor for each patient (higher = higher risk)
lp_train <- as.numeric(
  predict(cv_cox, newx = X, s = "lambda.min", type = "link")
)


library(dplyr)
# 2) Store time, status, and score
time   <- survtime
status <- outcome

risk_df <- data.frame(
  Patient_ID = data$Patient_ID[!indx],  # adjust to your row filter if needed
  time       = time,
  status     = status,
  risk_score = lp_train
)

# 3) Compute tertile cutoffs and assign groups
cutpoints <- quantile(risk_df$risk_score,
                      probs = c(0, 1/3, 2/3, 1),
                      na.rm = TRUE)

cutpoints

risk_df$risk_group3 <- cut(
  risk_df$risk_score,
  breaks = cutpoints,
  include.lowest = TRUE,
  labels = c("Low", "Intermediate", "High")
)

# Quick check
table(risk_df$risk_group3)

# 4) Save database with patient-level group membership
write.csv(
  risk_df[, c("Patient_ID", "time", "status", "risk_score", "risk_group3")],
  file = "CMML_Cox_RiskTertiles_PatientGroups.csv",
  row.names = FALSE
)


# coefs_df already has Variable and Coefficient, filtered to non-zero
# Add integer points scaled so largest |coef| ≈ 5 points
max_abs_coef <- max(abs(coefs_df$Coefficient))
coefs_df$Points <- round(5 * coefs_df$Coefficient / max_abs_coef)

# Save points table
write.csv(
  coefs_df[, c("Variable", "Coefficient", "HazardRatio", "Points")],
  file = "CMML_Cox_PointsTable.csv",
  row.names = FALSE
)

# Compute total points per patient
X_points <- X[, coefs_df$Variable, drop = FALSE]
total_points <- as.vector(X_points %*% coefs_df$Points)

risk_df$total_points <- total_points

# Tertiles based on total points (almost same grouping as risk_score)
pt_cutpoints <- quantile(total_points,
                         probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

risk_df$risk_group3_points <- cut(
  risk_df$total_points,
  breaks = pt_cutpoints,
  include.lowest = TRUE,
  labels = c("Low", "Intermediate", "High")
)

write.csv(
  risk_df[, c("Patient_ID", "time", "status",
              "risk_score", "total_points",
              "risk_group3", "risk_group3_points")],
  file = "CMML_Cox_RiskGroups_ScoreAndPoints.csv",
  row.names = FALSE
)






# rcorr.cens expects larger x = longer survival, so we use -lp
cobj <- rcorr.cens(-lp_train, y_surv)
apparent_cindex <- as.numeric(cobj["C Index"])

cat("\n==== APPARENT DISCRIMINATION ====\n")
cat("Harrell's C-index (apparent):", round(apparent_cindex, 3), "\n")



library(dplyr)
library(survminer)
library(ggplot2)

# Reuse your publication theme if already defined, or define one:
theme_publication <- function(base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.border     = element_rect(fill = NA, colour = "black", size = 0.5),
      panel.grid.major = element_line(colour = "grey90", size = 0.25),
      panel.grid.minor = element_blank(),
      axis.text        = element_text(size = rel(0.9), colour = "black"),
      axis.title       = element_text(size = rel(1.0), face = "bold"),
      axis.ticks       = element_line(colour = "black", size = 0.5),
      axis.line        = element_line(colour = "black", size = 0.5),
      legend.key       = element_rect(fill = "white", colour = NA),
      legend.title     = element_text(size = rel(0.9), face = "bold"),
      legend.text      = element_text(size = rel(0.8)),
      strip.background = element_rect(fill = "grey95", colour = "black", size = 0.5),
      strip.text       = element_text(size = rel(0.9), face = "bold"),
      plot.title       = element_text(size = rel(1.2), face = "bold",
                                      hjust = 0.5, margin = margin(b = 10)),
      plot.margin      = unit(c(1, 1, 1, 1), "lines")
    )
}


time = survtime
status = outcome

risk_df <- data.frame(
  time      = time,
  status    = status,
  risk_score = lp_train
)

# Tertiles of risk score for 3-group KM plot
risk_df$risk_group3 <- cut(
  risk_df$risk_score,
  breaks = quantile(risk_df$risk_score, probs = c(0, 1/3, 2/3, 1)),
  include.lowest = TRUE,
  labels = c("Low", "Intermediate", "High")
)

fit_km3 <- survfit(Surv(time, status) ~ risk_group3, data = risk_df)

fig_km3 <- ggsurvplot(
  fit_km3,
  data           = risk_df,
  risk.table     = TRUE,
  pval           = TRUE,
  conf.int       = TRUE,
  xlab           = "Follow-up time",
  ylab           = "Event-free survival",
  legend.title   = "Risk group",
  legend.labs    = levels(risk_df$risk_group3),
  palette        = c("#4575B4", "#FDB863", "#D73027"),
  ggtheme        = theme_publication()
)

fig_km3

ggsave("Figure_KM_RiskTertiles.pdf", fig_km3$plot,
       width = 8, height = 6, dpi = 300)


# fit_km3 already defined:
# fit_km3 <- survfit(Surv(time, status) ~ risk_group3, data = risk_df)

fig_km3_cum <- ggsurvplot(
  fit_km3,
  data         = risk_df,
  fun          = "event",
  risk.table   = TRUE,
  pval         = TRUE,
  pval.coord   = c(0.02, 0.90),   # (x, y) in data coordinates
  pval.size    = 4,
  conf.int     = TRUE,
  xlab         = "Follow-up time",
  ylab         = "Cumulative incidence of CMML progression",
  legend.title = "Risk group",
  legend.labs  = levels(risk_df$risk_group3),
  palette      = c("#4575B4", "#FDB863", "#D73027"),
  ggtheme      = theme_publication()
)

fig_km3_cum
# Save publication-quality versions
ggsave("Figure_KM_RiskTertiles_CumInc.pdf", fig_km3_cum$plot,
       width = 8, height = 6, dpi = 300)
ggsave("Figure_KM_RiskTertiles_CumInc.tiff", fig_km3_cum$plot,
       width = 8, height = 6, dpi = 300, compression = "lzw")



# Install once if needed:
# install.packages("timeROC")

library(timeROC)

lp_train <- as.numeric(
  predict(cv_cox, newx = X, s = "lambda.min", type = "link")
)
time   <- survtime
status <- outcome
# ==========================================
# TIME-DEPENDENT ROC AND AUC(t) FOR COX MODEL
# ==========================================

# Choose evaluation times (e.g., quartiles of event times)
event_times <- time[status == 1]
eval_times  <- quantile(event_times, probs = c(0.1,0.2,0.3,0.4, 0.5, 0.6, 0.7,0.8), na.rm = TRUE)
eval_times  <- as.numeric(eval_times)

cat("\nEvaluation times for AUC(t):", paste(round(eval_times, 2), collapse = ", "), "\n")

# Compute time-dependent ROC / AUC(t)
time_roc <- timeROC(
  T        = time,
  delta    = status,
  marker   = lp_train,   # linear predictor as risk score
  cause    = 1,          # event of interest
  weighting = "marginal",
  times    = eval_times,
  iid      = TRUE
)

# Extract AUC(t)
auc_df <- data.frame(
  time = time_roc$times,
  AUC  = time_roc$AUC
)

print(auc_df)

# --------------------------------------------------
# PUBLICATION-QUALITY PLOT: AUC(t) vs TIME
# --------------------------------------------------

fig_auc_time <- ggplot(auc_df, aes(x = time, y = AUC)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0.5, 1),
                     breaks = seq(0.5, 1.0, by = 0.1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "Time-dependent AUC for penalized Cox model",
    x     = "Time since baseline",
    y     = "AUC(t)"
  ) +
  theme_publication()

# Save as high-quality PDF and TIFF
ggsave("Figure_TimeDependent_AUC.pdf", fig_auc_time,
       width = 7, height = 5, dpi = 300)


# --------------------------------------------------
# OPTIONAL: ROC CURVE AT A SPECIFIC TIME POINT
# --------------------------------------------------

# Choose index for the evaluation time (e.g., median event time)
k <- 2  # 1 = 25th percentile, 2 = 50th, 3 = 75th

roc_df <- data.frame(
  FPR = time_roc$FP[, k],   # false positive rate
  TPR = time_roc$TP[, k]    # true positive rate
)
auc_at_t <- time_roc$AUC[k]

fig_roc_t <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", colour = "grey50") +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = paste0("Time-dependent ROC at t = ",
                   round(eval_times[k], 2)),
    x     = "False positive rate",
    y     = "True positive rate"
  ) +
  # Add C-statistics (and optionally AUC(t)) on the figure
  annotate(
    "text",
    x = 0.65, y = 0.15, hjust = 0,
    label = paste0(
      "C-index = ", sprintf("%.4f", apparent_cindex),
      "\nAUC(t)  = ", sprintf("%.4f", auc_at_t)
    ),
    size = 4
  ) +
  theme_publication()

ggsave("Figure_ROC_Timepoint.pdf", fig_roc_t,
       width = 7, height = 7, dpi = 300)
ggsave("Figure_ROC_Timepoint.tiff", fig_roc_t,
       width = 7, height = 7, dpi = 300, compression = "lzw")


