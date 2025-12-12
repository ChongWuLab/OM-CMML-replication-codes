# CMML Replication Codes 

This repository contains R scripts used for our manuscript. The goal of this README is to help readers reproduce the analyses and regenerate the key tables and figures.

> **Note:** Patient-level source data (CSV/Excel files) are not included in this repository. You will need access to the original data files described below and place them in the same directory as these scripts (or update paths accordingly).

---

## Repository Contents

- `comprehensive-cmml-analysis.r`  
  End-to-end CMML classification pipeline: feature engineering, unsupervised clustering, t‑SNE visualization, and elastic net classification with bootstrap inference and risk prediction for OM‑CMML patients.

- `CMLL_elasticnet_model_newAug28.r`  
  Binary (yes/no) CMML progression prediction using elastic net logistic regression, with extensive internal validation (bootstrap, optimism correction, repeated cross‑validation) and a points‑based clinical prediction model.

- `CMLL_elasticnet_model_newAug28-time-to-event.r`  
  Time‑to‑event (survival) modeling of CMML progression using penalized Cox regression (elastic net) and Kaplan–Meier curves by risk tertiles.

- `differential_gene_expression.r`  
  Transcriptomic pipeline for unsupervised clustering and differential expression (DESeq2), followed by pathway enrichment and publication‑quality figures (t‑SNE, volcano plots, MA plots, enrichment summaries).

- `CMML_Cox_Clinical_Prediction_Model_Variables.csv`  
  Example output: variables selected in the Cox elastic‑net time‑to‑event model.

- `Figure_KM_RiskTertiles.pdf`  
  Example output: Kaplan–Meier curves for risk tertiles from the time‑to‑event model.

---

## Software Requirements

- R (≥ 4.1 recommended)
- Suggested R packages (superset across scripts):
  - Data manipulation / plotting: `tidyverse`, `dplyr`, `ggplot2`, `cowplot`, `gridExtra`, `reshape2`, `ggrepel`, `patchwork`, `RColorBrewer`, `viridis`, `ggthemes`, `scales`
  - Machine learning / modeling: `caret`, `glmnet`, `randomForest`, `e1071`, `nnet`, `xgboost`, `kernlab`, `cluster`, `mclust`, `dbscan`, `vegan`, `Rtsne`, `ncvreg`, `hdm`, `hdi`, `selectiveInference`, `boot`, `PRROC`, `pROC`, `OptimalCutpoints`, `survival`, `Hmisc`, `survminer`
  - Transcriptomics / enrichment: `readxl`, `limma`, `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`, `pheatmap`, `ComplexHeatmap`, `circlize`

The scripts `comprehensive-cmml-analysis.r` and `differential_gene_expression.r` include helper code to install any missing packages automatically. For more controlled environments, you may prefer to install the listed packages manually via `install.packages()` and `BiocManager::install()`.

