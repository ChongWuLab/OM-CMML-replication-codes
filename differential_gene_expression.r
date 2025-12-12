#!/usr/bin/env Rscript

################################################################################
### Publication-Quality Analysis Script for MDS/CMML Transcriptomics
### Enhanced version with bug fixes and improved robustness
################################################################################

# --- 1. SETUP: LOAD LIBRARIES WITH INSTALLATION CHECK ---
required_packages <- c(
  "readxl", "dplyr", "ggplot2", "Rtsne", "ggrepel", "patchwork",
  "RColorBrewer", "viridis", "limma", "clusterProfiler", "org.Hs.eg.db",
  "ggthemes", "DESeq2", "pheatmap", "tidyr",
  "scales", "cowplot", "ComplexHeatmap", "circlize"
)

# Function to check and install packages
check_and_install <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "DESeq2", "ComplexHeatmap")) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }
}

# Load all required packages
check_and_install(required_packages)

# --- 2. PUBLICATION-QUALITY THEME AND COLOR SCHEMES ---

# Enhanced publication theme
theme_publication <- function(base_size = 14, base_family = "sans") {
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      # Title and text
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5,
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = rel(0.9), hjust = 0.5,
                                   margin = margin(b = 5)),
      text = element_text(color = "black"),
      
      # Background
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      
      # Panel
      panel.border = element_rect(fill = NA, colour = "black", size = 1),
      panel.grid.major = element_line(colour = "grey90", size = 0.25),
      panel.grid.minor = element_blank(),
      
      # Axes
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(size = rel(0.8), color = "black"),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", size = 0.5),
      
      # Legend
      legend.key = element_rect(fill = "white", colour = NA),
      legend.key.size = unit(0.6, "cm"),
      legend.spacing = unit(0.2, "cm"),
      legend.title = element_text(face = "bold", size = rel(0.9)),
      legend.text = element_text(size = rel(0.8)),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.box.background = element_rect(fill = "white", colour = "black", size = 0.5),
      
      # Facets
      strip.background = element_rect(fill = "grey95", colour = "black", size = 0.5),
      strip.text = element_text(face = "bold", size = rel(0.9)),
      
      # Margins
      plot.margin = unit(c(10, 10, 10, 10), "mm")
    )
}

# Set global theme
theme_set(theme_publication())

# Color palettes
# Colorblind-friendly qualitative palette
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Significance colors for volcano plots
sig_colors <- c("Upregulated" = "#D55E00",
                "Downregulated" = "#0072B2",
                "Not Significant" = "grey70")

# --- 3. HELPER FUNCTIONS ---

# Function to inspect Excel file structure
inspect_excel_structure <- function(filename) {
  cat("\n=== Excel File Structure Inspection ===\n")
  
  # Read first 20 rows to understand structure
  preview <- read_excel(filename, n_max = 20, col_names = FALSE)
  
  cat("\nFirst column (should contain metadata field names):\n")
  first_col <- as.character(unlist(preview[, 1]))
  for (i in 1:min(10, length(first_col))) {
    if (!is.na(first_col[i]) && first_col[i] != "") {
      cat("  Row", i, ":", first_col[i], "\n")
    }
  }
  
  cat("\nChecking row 1 for sample IDs:\n")
  first_row <- as.character(unlist(preview[1, ]))
  # Find where actual sample IDs start
  sample_start <- which(!is.na(first_row) & first_row != "" &
                       !first_row %in% c("Sample_ID", "A", "B", "C", letters, LETTERS))[1]
  
  if (!is.na(sample_start)) {
    sample_ids <- first_row[sample_start:length(first_row)]
    sample_ids <- sample_ids[!is.na(sample_ids) & sample_ids != ""]
    cat("  Sample IDs start at column", sample_start, "\n")
    cat("  First few samples:", paste(head(sample_ids, 5), collapse = ", "), "...\n")
    cat("  Total samples detected:", length(sample_ids), "\n")
  }
  
  # Check for RNA-seq data row
  rna_row <- which(apply(preview, 1, function(x) {
    any(grepl("RNA.seq|RNA-seq data", x, ignore.case = TRUE))
  }))[1]
  
  if (!is.na(rna_row)) {
    cat("\n'RNA-seq data' header found at row:", rna_row, "\n")
    cat("Expression data should start at row:", rna_row + 1, "\n")
  }
  
  # Check data types in a few cells
  if (!is.na(sample_start) && nrow(preview) > 10) {
    cat("\nChecking data types in expression matrix area:\n")
    test_val <- preview[11, sample_start]
    cat("  Sample value at row 11:", test_val, "(type:", class(test_val), ")\n")
  }
  
  cat("\n")
}


# Function to filter and transform expression data
prepare_expression_data <- function(expr_matrix, n_top_var = 5000, min_expr = 1) {
  # Remove genes with low expression
  keep <- rowSums(expr_matrix >= min_expr) >= 3
  expr_filtered <- expr_matrix[keep, ]
  cat("Genes after low expression filter:", nrow(expr_filtered), "\n")
  
  # Log2 transformation
  expr_log <- log2(expr_filtered + 1)
  
  # Select top variable genes
  gene_vars <- apply(expr_log, 1, var, na.rm = TRUE)
  top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(n_top_var, length(gene_vars))])
  
  cat("Selected top", length(top_var_genes), "variable genes for clustering\n")
  
  return(list(
    full = expr_log,
    variable = expr_log[top_var_genes, ]
  ))
}



# --- 4. MAIN ANALYSIS ---

input_file = "Copy of merged_data_Selected by GMB 2025.xlsx"
# Get filename from user or use default
if (exists("input_file") && !is.null(input_file)) {
  filename <- input_file
} else {
  # List Excel files in current directory
  excel_files <- list.files(pattern = "\\.(xlsx|xls)$", ignore.case = TRUE)
  if (length(excel_files) > 0) {
    cat("Excel files found in current directory:\n")
    for (i in 1:length(excel_files)) {
      cat(i, ":", excel_files[i], "\n")
    }
    cat("\nUsing first file. To use a different file, set 'input_file' variable.\n")
    filename <- excel_files[1]
  } else {
    stop("No Excel files found. Please set 'input_file' variable to your data file path.")
  }
}

# Alternative function to read Excel with automatic format detection
read_excel_auto <- function(filename) {
  cat("Attempting to read Excel file with automatic format detection...\n")
  
  # Read all data
  all_data <- read_excel(filename, col_names = FALSE)
  
  # Find where RNA-seq data starts (look for "RNA-seq data" or similar)
  rna_seq_row <- which(apply(all_data, 1, function(x) {
    any(grepl("RNA.seq|gene.*type|gene.*name|ENSG", x, ignore.case = TRUE))
  }))[1]
  
  if (is.na(rna_seq_row)) {
    # If not found, assume it starts after 8 rows of metadata
    rna_seq_row <- 9
  }
  
  cat("Detected RNA-seq data starting at row:", rna_seq_row, "\n")
  
  # Extract metadata
  metadata_rows <- rna_seq_row - 1
  metadata_raw <- all_data[1:metadata_rows, ]
  
  # Get sample IDs from first row
  sample_row <- as.character(unlist(metadata_raw[1, ]))
  sample_start <- which(!is.na(sample_row) & sample_row != "Sample_ID" & sample_row != "")[1]
  sample_ids <- sample_row[sample_start:length(sample_row)]
  sample_ids <- sample_ids[!is.na(sample_ids) & sample_ids != ""]
  
  cat("Found", length(sample_ids), "samples\n")
  
  # Create metadata dataframe
  metadata <- data.frame(Sample_ID = sample_ids)
  
  # Add other metadata fields
  metadata_names <- as.data.frame(metadata_raw[1:metadata_rows, 1])[,1]
  
  for (i in 2:metadata_rows) {
    field_name <- metadata_names[i]
    if (!is.na(field_name) && field_name != "") {
      field_values <- as.character(unlist(metadata_raw[i, sample_start:(sample_start + length(sample_ids) - 1)]))
      metadata[[field_name]] <- field_values
    }
  }
  
  # Clean column names
  colnames(metadata) <- gsub("[^[:alnum:]_]", "_", colnames(metadata))
  colnames(metadata) <- gsub("biTET2_SRSF2", "biTET2_SRSF2", colnames(metadata))
  rownames(metadata) <- metadata$Sample_ID
  
  # Extract expression data
  # First, find where actual numeric data starts
  expr_start <- rna_seq_row + 1
  
  # Read expression matrix
  expr_data <- all_data[expr_start:nrow(all_data), sample_start:(sample_start + length(sample_ids) - 1)]
  expression_matrix <- as.matrix(expr_data)
  
  # Get gene names
  gene_col <- all_data[expr_start:nrow(all_data), 1]
  gene_names <- as.character(unlist(gene_col))
  
  # Clean gene names
  gene_names <- gene_names[!is.na(gene_names) & gene_names != ""]
  
  if (length(gene_names) == nrow(expression_matrix)) {
    rownames(expression_matrix) <- gene_names
  } else {
    warning("Gene names don't match expression matrix rows. Using generic names.")
    rownames(expression_matrix) <- paste0("Gene_", 1:nrow(expression_matrix))
  }
  
  colnames(expression_matrix) <- sample_ids
  
  # Convert to numeric
  class(expression_matrix) <- "numeric"
  
  # Remove rows with all NAs
  keep_rows <- rowSums(!is.na(expression_matrix)) > 0
  expression_matrix <- expression_matrix[keep_rows, ]
  
  return(list(metadata = metadata, expression = expression_matrix))
}

# Main data loading section with fallback
data_list <- read_excel_auto(filename)

expression_matrix = data_list$expression

# Prepare expression data
expr_data <- prepare_expression_data(expression_matrix)
expr_full_log <- expr_data$full
expr_var_log <- expr_data$variable

metadata = data_list$metadata
metadata[grepl("ctrl",metadata$Patient_ID),"Genomic_group"] = "Control"
metadata[grepl("ctrl",metadata$Patient_ID),"Genomic_cluster"] = "Control"
################################################################################
### PART 1: PUBLICATION-QUALITY CLUSTERING ANALYSIS (REVISED)
################################################################################

cat("\n--- PART 1: Clustering Analysis ---\n")

# Create output directory
dir.create("Publication_Figures", showWarnings = FALSE)

# --- 1A. Hierarchical Clustering ---
cat("Performing hierarchical clustering...\n")

# Prepare data for clustering (samples as columns, genes as rows)
# Using the variable genes for clustering
clustering_data <- expr_var_log

# Calculate distance matrix between samples
sample_dist <- dist(t(clustering_data), method = "euclidean")

# Perform hierarchical clustering
hc <- hclust(sample_dist, method = "ward.D2")

# Plot dendrogram to help determine optimal number of clusters
pdf("Publication_Figures/Dendrogram_Samples.pdf", width = 12, height = 8)
plot(hc, main = "Hierarchical Clustering of Samples",
     xlab = "Samples", ylab = "Distance",
     cex = 0.8, hang = -1)
dev.off()

# Determine optimal number of clusters using elbow method
# Calculate within-cluster sum of squares for different k values
wss <- sapply(2:10, function(k) {
  clusters <- cutree(hc, k = k)
  sum(sapply(unique(clusters), function(i) {
    samples_in_cluster <- names(clusters)[clusters == i]
    if (length(samples_in_cluster) > 1) {
      cluster_data <- clustering_data[, samples_in_cluster]
      sum(dist(t(cluster_data))^2) / (2 * ncol(cluster_data))
    } else {
      0
    }
  }))
})

# Plot elbow curve
pdf("Publication_Figures/Elbow_Plot_Clustering.pdf", width = 8, height = 6)
plot(2:10, wss, type = "b", pch = 19,
     xlab = "Number of Clusters",
     ylab = "Within-Cluster Sum of Squares",
     main = "Elbow Method for Optimal Clusters")
dev.off()

# Based on elbow method and typical expectations, let's try k=3 to k=5
# We'll visualize different k values and see which makes most biological sense
k_values <- c(3, 4, 5)

# --- 1B. t-SNE Visualization of Clusters ---
cat("Running t-SNE analysis for cluster visualization...\n")

# Prepare data for t-SNE (samples as rows)
data_for_tsne <- t(clustering_data)

# Calculate appropriate perplexity
n_samples <- nrow(data_for_tsne)
perplexity <- min(30, floor((n_samples - 1) / 3))
cat("Using perplexity:", perplexity, "\n")

# Run t-SNE
set.seed(42)
tsne_result <- tryCatch({
  Rtsne(data_for_tsne,
        dims = 2,
        perplexity = perplexity,
        check_duplicates = FALSE,
        pca = TRUE,
        pca_center = TRUE,
        pca_scale = FALSE,
        verbose = TRUE,
        max_iter = 1000)
}, error = function(e) {
  cat("Error in t-SNE, trying with different parameters...\n")
  Rtsne(data_for_tsne,
        dims = 2,
        perplexity = 10,
        check_duplicates = FALSE,
        pca = FALSE,
        verbose = TRUE)
})

# Create base plotting dataframe
tsne_df <- data.frame(
  tSNE1 = tsne_result$Y[,1],
  tSNE2 = tsne_result$Y[,2],
  Sample_ID = rownames(data_for_tsne)
)

# Add hierarchical clustering results for different k values
for (k in k_values) {
  cluster_col <- paste0("HC_k", k)
  tsne_df[[cluster_col]] <- factor(cutree(hc, k = k)[tsne_df$Sample_ID])
}

# Merge with metadata
tsne_df <- merge(tsne_df, metadata, by = "Sample_ID")

# --- 1C. Create Multi-Panel Figure Showing Clusters and Metadata ---

# Helper function for t-SNE plots with dual coloring
create_tsne_dual_plot <- function(data, cluster_var, metadata_var, title, k_value) {
  # Create color mapping for clusters
  n_clusters <- length(unique(data[[cluster_var]]))
  cluster_colors <- brewer.pal(max(3, n_clusters), "Set1")[1:n_clusters]
  
  # Create shape mapping for metadata
  metadata_values <- unique(data[[metadata_var]])
  n_metadata <- length(metadata_values)
  shape_values <- c(21, 22, 23, 24, 25)[1:min(5, n_metadata)]
  
  p <- ggplot(data, aes_string(x = "tSNE1", y = "tSNE2")) +
    geom_point(aes_string(fill = cluster_var, shape = paste0("`", metadata_var, "`")),
               size = 6, color = "black", stroke = 0.8, alpha = 0.8) +
    scale_fill_manual(values = cluster_colors, name = paste("Cluster (k =", k_value, ")")) +
    scale_shape_manual(values = shape_values, name = gsub("_", " ", metadata_var)) +
    labs(x = "t-SNE Dimension 1",
         y = "t-SNE Dimension 2",
         title = title) +
    theme_publication() +
    theme(legend.position = "right",
          legend.box = "vertical") +
    guides(fill = guide_legend(order = 1, override.aes = list(shape = 21)),
           shape = guide_legend(order = 2, override.aes = list(size = 4)))
  
  return(p)
}

# Create individual plots for each k value and metadata combination
plot_list <- list()
plot_index <- 1

# Check which metadata columns are available
metadata_cols <- c("Diagnostic_group", "Genomic_cluster", "Genomic_group")
available_metadata <- metadata_cols[metadata_cols %in% colnames(tsne_df)]

if (length(available_metadata) == 0) {
  # Try to find similar column names
  cat("Looking for metadata variables...\n")
  possible_vars <- grep("Diagnostic|Genomic|cluster|group", colnames(tsne_df),
                        value = TRUE, ignore.case = TRUE)
  cat("Found variables:", paste(possible_vars, collapse = ", "), "\n")
  # Filter out the HC_k columns we just added
  available_metadata <- possible_vars[!grepl("^HC_k", possible_vars)][1:min(3, length(possible_vars))]
}

# Create plots for k=4 (middle value) with each metadata variable
optimal_k <- 5
cluster_col <- paste0("HC_k", optimal_k)

for (meta_var in available_metadata) {
  plot_title <- paste("Hierarchical Clustering vs", gsub("_", " ", meta_var))
  
  plot_list[[plot_index]] <- create_tsne_dual_plot(
    tsne_df, cluster_col, meta_var, plot_title, optimal_k
  )
  plot_index <- plot_index + 1
}

# --- 1D. Create Comparison Plots ---

# Function to calculate cluster purity with respect to metadata groups
calculate_cluster_purity <- function(clusters, metadata_groups) {
  purity <- 0
  n <- length(clusters)
  
  for (i in unique(clusters)) {
    cluster_samples <- which(clusters == i)
    if (length(cluster_samples) > 0) {
      # Find most common metadata group in this cluster
      group_counts <- table(metadata_groups[cluster_samples])
      max_count <- max(group_counts)
      purity <- purity + max_count
    }
  }
  
  return(purity / n)
}

# Calculate purity for each k and metadata combination
purity_results <- data.frame()

for (k in k_values) {
  cluster_col <- paste0("HC_k", k)
  clusters <- tsne_df[[cluster_col]]
  
  for (meta_var in available_metadata) {
    metadata_groups <- tsne_df[[meta_var]]
    purity <- calculate_cluster_purity(clusters, metadata_groups)
    
    purity_results <- rbind(purity_results, data.frame(
      k = k,
      metadata = meta_var,
      purity = purity
    ))
  }
}

# Create purity comparison plot
purity_plot <- ggplot(purity_results, aes(x = factor(k), y = purity,
                                           fill = metadata, group = metadata)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = cb_palette, name = "Metadata Variable") +
  labs(x = "Number of Clusters (k)",
       y = "Cluster Purity",
       title = "Cluster Purity Analysis",
       subtitle = "Agreement between hierarchical clusters and metadata groups") +
  theme_publication() +
  theme(legend.position = "bottom")

plot_list[[plot_index]] <- purity_plot

# --- 1E. Create Silhouette Analysis ---

# Calculate silhouette scores for different k values
library(cluster)

silhouette_plot_data <- data.frame()

for (k in k_values) {
  clusters <- cutree(hc, k = k)
  sil <- silhouette(clusters, sample_dist)
  avg_sil <- mean(sil[, 3])
  
  silhouette_plot_data <- rbind(silhouette_plot_data, data.frame(
    k = k,
    avg_silhouette = avg_sil
  ))
}

silhouette_plot <- ggplot(silhouette_plot_data, aes(x = factor(k), y = avg_silhouette)) +
  geom_bar(stat = "identity", fill = cb_palette[1], alpha = 0.8) +
  geom_text(aes(label = round(avg_silhouette, 3)), vjust = -0.5) +
  labs(x = "Number of Clusters (k)",
       y = "Average Silhouette Score",
       title = "Silhouette Analysis",
       subtitle = "Higher scores indicate better-defined clusters") +
  theme_publication() +
  ylim(0, max(silhouette_plot_data$avg_silhouette) * 1.1)

plot_list[[plot_index + 1]] <- silhouette_plot

# --- 1F. Create Sample Assignment Table ---

# Create a summary table showing cluster assignments
cluster_assignments <- data.frame(
  Sample_ID = tsne_df$Sample_ID,
  HC_Cluster_k4 = tsne_df$HC_k4
)

# Add metadata columns
for (meta_var in available_metadata) {
  cluster_assignments[[meta_var]] <- tsne_df[[meta_var]]
}

# Save cluster assignments
write.csv(cluster_assignments,
          "Publication_Figures/Cluster_Assignments.csv",
          row.names = FALSE)

# --- 1G. Combine Plots into Publication Figure ---

# Main figure with clustering results
figure1_main <- wrap_plots(plot_list[1:length(available_metadata)], ncol = 1) +
  plot_annotation(
    title = 'Transcriptomic Clustering of MDS and CMML Samples',
    subtitle = 'Hierarchical clustering (Ward.D2) based on top 5000 variable genes',
    tag_levels = 'A',
    theme = theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      plot.subtitle = element_text(size = 16, hjust = 0.5)
    )
  )

# Supplementary figure with analysis metrics
figure1_supp <- wrap_plots(
  list(purity_plot, silhouette_plot),
  ncol = 2
) +
  plot_annotation(
    title = 'Clustering Quality Metrics',
    tag_levels = 'A',
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
    )
  )

# Save figures
ggsave("Publication_Figures/Figure1_Clustering_Analysis.pdf",
       figure1_main, width = 10, height = 6 * length(available_metadata),
       device = cairo_pdf)
ggsave("Publication_Figures/Figure1_Clustering_Analysis.png",
       figure1_main, width = 10, height = 6 * length(available_metadata),
       dpi = 300)

ggsave("Publication_Figures/Figure1_Supp_Clustering_Metrics.pdf",
       figure1_supp, width = 12, height = 6, device = cairo_pdf)

# --- 1H. Create detailed cluster characterization ---


# Option 1: Create a detailed summary table
create_cluster_summary <- function(tsne_df, cluster_column, metadata_vars) {
  clusters <- tsne_df[[cluster_column]]
  summary_list <- list()
  
  for (i in sort(unique(clusters))) {
    cluster_samples <- tsne_df$Sample_ID[clusters == i]
    
    # Start with basic info
    row_data <- data.frame(
      Cluster = i,
      N_Samples = length(cluster_samples),
      Sample_IDs = paste(cluster_samples, collapse = "; ")
    )
    
    # Add metadata breakdowns
    for (meta_var in metadata_vars) {
      meta_table <- table(tsne_df[[meta_var]][clusters == i])
      
      # Create columns for each unique value in this metadata variable
      for (meta_value in names(meta_table)) {
        col_name <- paste0(meta_var, "_", meta_value)
        count <- meta_table[meta_value]
        percent <- round(100 * count / sum(meta_table), 1)
        row_data[[col_name]] <- paste0(count, " (", percent, "%)")
      }
    }
    
    summary_list[[i]] <- row_data
  }
  
  # Combine all rows
  summary_df <- do.call(rbind, summary_list)
  return(summary_df)
}

# Option 2: Create a long-format detailed table
create_cluster_details <- function(tsne_df, cluster_column, metadata_vars) {
  clusters <- tsne_df[[cluster_column]]
  details_list <- list()
  
  for (i in sort(unique(clusters))) {
    cluster_samples <- tsne_df$Sample_ID[clusters == i]
    n_samples <- length(cluster_samples)
    
    # Add a summary row
    details_list[[length(details_list) + 1]] <- data.frame(
      Cluster = i,
      Category = "Summary",
      Variable = "Total",
      Value = n_samples,
      Count = n_samples,
      Percentage = 100,
      Samples = paste(cluster_samples, collapse = "; ")
    )
    
    # Add metadata breakdowns
    for (meta_var in metadata_vars) {
      meta_table <- table(tsne_df[[meta_var]][clusters == i])
      
      for (meta_value in names(meta_table)) {
        count <- meta_table[meta_value]
        percent <- round(100 * count / sum(meta_table), 1)
        
        # Get sample IDs for this specific combination
        samples_in_group <- tsne_df$Sample_ID[clusters == i &
                                               tsne_df[[meta_var]] == meta_value]
        
        details_list[[length(details_list) + 1]] <- data.frame(
          Cluster = i,
          Category = gsub("_", " ", meta_var),
          Variable = meta_value,
          Value = paste0(count, " (", percent, "%)"),
          Count = count,
          Percentage = percent,
          Samples = paste(samples_in_group, collapse = "; ")
        )
      }
    }
  }
  
  details_df <- do.call(rbind, details_list)
  return(details_df)
}

# Apply both options to your data
optimal_clusters <- tsne_df$HC_k5

# Create summary table (wide format)
cluster_summary <- create_cluster_summary(tsne_df, "HC_k5", available_metadata)
write.csv(cluster_summary,
          "Publication_Figures/Cluster_Summary_k5.csv",
          row.names = FALSE)

# Create detailed table (long format)
cluster_details <- create_cluster_details(tsne_df, "HC_k5", available_metadata)
write.csv(cluster_details,
          "Publication_Figures/Cluster_Details_k5.csv",
          row.names = FALSE)

# Also create a simple cluster membership table
cluster_membership <- data.frame(
  Sample_ID = tsne_df$Sample_ID,
  Cluster = tsne_df$HC_k5
)

# Add all metadata columns
for (meta_var in available_metadata) {
  cluster_membership[[meta_var]] <- tsne_df[[meta_var]]
}

# Sort by cluster
cluster_membership <- cluster_membership[order(cluster_membership$Cluster), ]
write.csv(cluster_membership,
          "Publication_Figures/Cluster_Membership_k5.csv",
          row.names = FALSE)

# Print the same output as before but also save it
cat("\nCluster Characterization (k=5):\n")
cat("================================\n")

# Create a text file with the console output
sink("Publication_Figures/Cluster_Characterization_k5.txt")

for (i in sort(unique(optimal_clusters))) {
  cat("\nCluster", i, ":\n")
  cluster_samples <- tsne_df$Sample_ID[optimal_clusters == i]
  cat("  Number of samples:", length(cluster_samples), "\n")
  cat("  Sample IDs:", paste(cluster_samples, collapse = ", "), "\n")
  
  for (meta_var in available_metadata) {
    cat("  ", gsub("_", " ", meta_var), ":\n", sep = "")
    meta_table <- table(tsne_df[[meta_var]][optimal_clusters == i])
    for (j in 1:length(meta_table)) {
      cat("    ", names(meta_table)[j], ": ", meta_table[j],
          " (", round(100 * meta_table[j] / sum(meta_table), 1), "%)\n", sep = "")
    }
  }
}

sink()  # Close the text file

# Also print to console
for (i in sort(unique(optimal_clusters))) {
  cat("\nCluster", i, ":\n")
  cluster_samples <- tsne_df$Sample_ID[optimal_clusters == i]
  cat("  Number of samples:", length(cluster_samples), "\n")
  
  for (meta_var in available_metadata) {
    cat("  ", gsub("_", " ", meta_var), ":\n", sep = "")
    meta_table <- table(tsne_df[[meta_var]][optimal_clusters == i])
    for (j in 1:length(meta_table)) {
      cat("    ", names(meta_table)[j], ": ", meta_table[j],
          " (", round(100 * meta_table[j] / sum(meta_table), 1), "%)\n", sep = "")
    }
  }
}

cat("\n\nCluster characterization saved to multiple formats:\n")
cat("- Cluster_Summary_k5.csv (wide format summary)\n")
cat("- Cluster_Details_k5.csv (long format with all details)\n")
cat("- Cluster_Membership_k5.csv (sample-level assignments)\n")
cat("- Cluster_Characterization_k5.txt (formatted text output)\n")

cat("\nFigure 1 saved successfully\n")
cat("Cluster assignments saved to: Cluster_Assignments.csv\n")




################################################################################
### PART 2: DIFFERENTIAL EXPRESSION ANALYSIS WITH DESeq2 (FIXED)
################################################################################

cat("\n--- PART 2: Differential Expression Analysis ---\n")

# --- 2A. Prepare Data for DESeq2 ---

# First, we need the raw counts (not log-transformed)
# DESeq2 requires raw counts, not normalized or log-transformed data
raw_counts <- expression_matrix  # This should be the original count matrix

# Ensure counts are integers
raw_counts <- round(raw_counts)

# Find the biTET2/SRSF2 column
bitet_col <- grep("biTET2|TET2.*SRSF2|SRSF2", colnames(metadata),
                  ignore.case = TRUE, value = TRUE)[1]

if (is.na(bitet_col)) {
  possible_names <- c("biTET2/SRSF2", "biTET2_SRSF2", "biTET2SRSF2")
  for (name in possible_names) {
    if (name %in% colnames(metadata)) {
      bitet_col <- name
      break
    }
  }
}

diag_col <- grep("Diagnostic", colnames(metadata),
                 ignore.case = TRUE, value = TRUE)[1]

# Find Age and Sex columns
age_col <- grep("Age|age", colnames(metadata), ignore.case = TRUE, value = TRUE)[1]
sex_col <- grep("Sex|Gender|gender", colnames(metadata), ignore.case = TRUE, value = TRUE)[1]

cat("Using columns:\n")
cat("  Diagnostic column:", diag_col, "\n")
cat("  biTET2/SRSF2 column:", bitet_col, "\n")
cat("  Age column:", age_col, "\n")
cat("  Sex column:", sex_col, "\n\n")

# Create DE groups
metadata$DE_Group <- with(metadata, {
  bitet_status <- metadata[[bitet_col]]
  diag_group <- metadata[[diag_col]]
  
  case_when(
    grepl("Yes|TRUE|1", bitet_status, ignore.case = TRUE) &
      grepl("MDS", diag_group) ~ "MDS_biTET2_SRSF2",
    grepl("Yes|TRUE|1", bitet_status, ignore.case = TRUE) &
      grepl("CMML", diag_group) ~ "CMML_biTET2_SRSF2",
    grepl("No|FALSE|0", bitet_status, ignore.case = TRUE) &
      grepl("MDS", diag_group) ~ "MDS_Other",
    grepl("No|FALSE|0", bitet_status, ignore.case = TRUE) &
      grepl("CMML", diag_group) ~ "CMML_Other",
    TRUE ~ "Exclude"
  )
})

# Print group sizes
cat("\nDE Group sizes:\n")
print(table(metadata$DE_Group))

# Filter samples
keep_samples <- metadata$DE_Group != "Exclude"
metadata_de <- metadata[keep_samples, ]

# Handle duplicate gene names BEFORE filtering
if (any(duplicated(rownames(raw_counts)))) {
  cat("\nHandling duplicate gene names...\n")
  dup_genes <- rownames(raw_counts)[duplicated(rownames(raw_counts))]
  cat("Found", length(unique(dup_genes)), "duplicated gene names\n")
  
  # Option 1: Make unique by adding suffix
  rownames(raw_counts) <- make.unique(rownames(raw_counts), sep = "_")
}

# Now filter samples
raw_counts_de <- raw_counts[, keep_samples]

# Ensure metadata is in the same order as count matrix
metadata_de <- metadata_de[match(colnames(raw_counts_de), metadata_de$Sample_ID), ]
rownames(metadata_de) <- metadata_de$Sample_ID

# Clean up Age - center and scale to avoid collinearity issues
if (!is.na(age_col)) {
  age_values <- as.numeric(gsub("[^0-9.]", "", metadata_de[[age_col]]))
  
  # Check for missing ages
  if (any(is.na(age_values))) {
    cat("Warning: Found", sum(is.na(age_values)), "missing age values. Imputing with median.\n")
    median_age <- median(age_values, na.rm = TRUE)
    age_values[is.na(age_values)] <- median_age
  }
  
  # Center and scale age
  metadata_de$Age_scaled <- scale(age_values, center = TRUE, scale = TRUE)[,1]
  cat("Age statistics (before scaling):\n")
  cat("  Mean:", round(mean(age_values), 1), "\n")
  cat("  SD:", round(sd(age_values), 1), "\n")
  cat("  Range:", round(min(age_values), 1), "-", round(max(age_values), 1), "\n")
} else {
  cat("Warning: Age column not found. Creating dummy age variable.\n")
  metadata_de$Age_scaled <- 0  # Centered dummy variable
}

# Clean up Sex and convert to factor with specified levels
if (!is.na(sex_col)) {
  sex_values <- toupper(substr(metadata_de[[sex_col]], 1, 1))
  sex_values[!sex_values %in% c("M", "F")] <- "Unknown"
  
  # Convert to factor with explicit levels
  metadata_de$Sex_factor <- factor(sex_values, levels = c("F", "M", "Unknown"))
  
  cat("Sex distribution:\n")
  print(table(metadata_de$Sex_factor))
} else {
  cat("Warning: Sex column not found. Creating dummy sex variable.\n")
  metadata_de$Sex_factor <- factor("Unknown", levels = c("F", "M", "Unknown"))
}

# Convert DE_Group to factor with explicit levels
# Set reference level (the one that will be the baseline in comparisons)
de_groups <- unique(metadata_de$DE_Group)
de_groups <- de_groups[de_groups != "Exclude"]

# Set MDS_Other as reference (or choose your preferred reference)
if ("MDS_Other" %in% de_groups) {
  metadata_de$DE_Group_factor <- factor(metadata_de$DE_Group,
                                        levels = c("MDS_Other",
                                                  setdiff(de_groups, "MDS_Other")))
} else {
  metadata_de$DE_Group_factor <- factor(metadata_de$DE_Group)
}

cat("\nDE Group factor levels:\n")
print(levels(metadata_de$DE_Group_factor))

# --- 2B. DESeq2 Analysis ---
cat("\nPerforming differential expression analysis with DESeq2...\n")

# Remove any genes with all zero counts
non_zero_genes <- rowSums(raw_counts_de) > 0
raw_counts_de <- raw_counts_de[non_zero_genes, ]
cat("Removed", sum(!non_zero_genes), "genes with zero counts across all samples\n")

# Create DESeq2 dataset with properly formatted design
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts_de,
  colData = metadata_de,
  design = ~ Age_scaled + Sex_factor + DE_Group_factor
)

# Pre-filtering: remove genes with very low counts
keep_genes <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep_genes,]
cat("Genes after filtering:", nrow(dds), "\n")

# Check design matrix for any issues
cat("\nDesign matrix preview (first 5 samples):\n")
design_matrix <- model.matrix(~ Age_scaled + Sex_factor + DE_Group_factor, data = metadata_de)
print(head(design_matrix, 5))

# Run DESeq2
dds <- DESeq(dds)

# Extract normalized counts for later use
normalized_counts <- counts(dds, normalized = TRUE)

# --- 2C. Define Contrasts and Get Results ---

# Get the available result names to verify our contrasts
cat("\nAvailable contrasts in the model:\n")
print(resultsNames(dds))

# Update contrast definitions to use the factor column names
# The correct format for contrasts is: c("variable_name", "level1", "level2")
contrasts_list <- list(
  "MDS_biTET2_vs_CMML_biTET2" = c("DE_Group_factor", "MDS_biTET2_SRSF2", "CMML_biTET2_SRSF2"),
  "MDS_biTET2_vs_MDS_Other" = c("DE_Group_factor", "MDS_biTET2_SRSF2", "MDS_Other"),
  "CMML_biTET2_vs_CMML_Other" = c("DE_Group_factor", "CMML_biTET2_SRSF2", "CMML_Other")
)

# Set thresholds
FDR_THRESHOLD <- 0.05
LOGFC_THRESHOLD <- 1.0

# Process each comparison
volcano_plots <- list()
de_results_list <- list()

for (comp_name in names(contrasts_list)) {
  contrast <- contrasts_list[[comp_name]]
  cat("\nProcessing:", comp_name, "\n")
  
  # Get results
  res <- results(dds, contrast = contrast, alpha = FDR_THRESHOLD)
  
  # Convert to data frame and add gene names
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  
  # Add significance categories
  res_df$Significance <- case_when(
    res_df$padj < FDR_THRESHOLD & res_df$log2FoldChange > LOGFC_THRESHOLD ~ "Upregulated",
    res_df$padj < FDR_THRESHOLD & res_df$log2FoldChange < -LOGFC_THRESHOLD ~ "Downregulated",
    TRUE ~ "Not Significant"
  )
  
  # Remove NA p-values
  res_df <- res_df[!is.na(res_df$padj), ]
  
  # Store results
  de_results_list[[comp_name]] <- res_df
  
  # Count significant genes
  n_up <- sum(res_df$Significance == "Upregulated", na.rm = TRUE)
  n_down <- sum(res_df$Significance == "Downregulated", na.rm = TRUE)
  n_total <- n_up + n_down
  
  cat("  Upregulated:", n_up, "\n")
  cat("  Downregulated:", n_down, "\n")
  cat("  Total DE:", n_total, "\n")
  
  # Select top genes to label (by p-value and fold change)
  top_genes <- res_df %>%
    filter(Significance != "Not Significant") %>%
    mutate(score = -log10(padj) * abs(log2FoldChange)) %>%
    arrange(desc(score)) %>%
    head(20)
  
  # If we have ENSG IDs, try to clean them up
  if (any(grepl("ENSG", res_df$Gene))) {
    # Remove version numbers from ENSG IDs for cleaner labels
    top_genes$Gene_label <- gsub("\\.\\d+$", "", top_genes$Gene)
  } else {
    top_genes$Gene_label <- top_genes$Gene
  }
  
  # Create enhanced volcano plot
  volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    # Add background grid
    theme_publication() +
    
    # Add points with different colors
    geom_point(data = subset(res_df, Significance == "Not Significant"),
               aes(color = Significance), size = 1, alpha = 0.4) +
    geom_point(data = subset(res_df, Significance != "Not Significant"),
               aes(color = Significance), size = 2, alpha = 0.8) +
    
    # Color scale
    scale_color_manual(values = sig_colors, name = NULL) +
    
    # Add threshold lines
    geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
               linetype = "dashed", color = "grey40", size = 0.5, alpha = 0.8) +
    geom_hline(yintercept = -log10(FDR_THRESHOLD),
               linetype = "dashed", color = "grey40", size = 0.5, alpha = 0.8) +
    
    # Add gene labels with improved positioning
    geom_text_repel(data = top_genes,
                    aes(label = Gene_label),
                    size = 3,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    segment.color = 'grey50',
                    segment.size = 0.3,
                    segment.alpha = 0.8,
                    max.overlaps = 30,
                    force = 2,
                    min.segment.length = 0) +
    
    # Labels
    labs(title = gsub("_", " ", comp_name),
         subtitle = bquote(italic(n)[DE] ~ "=" ~ .(n_total) ~
                          " (" * uparrow * .(n_up) ~
                          ", " * downarrow * .(n_down) * ")"),
         x = expression(Log[2] ~ "Fold Change"),
         y = expression(-Log[10] ~ italic(P)[adjusted])) +
    
    # Theme adjustments
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box.background = element_blank(),
          panel.grid.major = element_line(colour = "grey90", size = 0.25),
          panel.grid.minor = element_blank()) +
    
    # Set axis limits to ensure symmetry
    coord_cartesian(xlim = c(-max(abs(res_df$log2FoldChange), na.rm = TRUE),
                             max(abs(res_df$log2FoldChange), na.rm = TRUE)))
  
  volcano_plots[[comp_name]] <- volcano
  
  # Save individual volcano plot
  ggsave(paste0("Publication_Figures/Volcano_", comp_name, ".pdf"),
         volcano, width = 8, height = 8, device = cairo_pdf)
  ggsave(paste0("Publication_Figures/Volcano_", comp_name, ".png"),
         volcano, width = 8, height = 8, dpi = 300, bg = "white")
  
  # Save DE results with additional statistics
  # FIX: Use dplyr::select to avoid conflicts
  results_output <- res_df %>%
    arrange(padj) %>%
    dplyr::select(Gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, Significance)
  
  write.csv(results_output,
            paste0("Publication_Figures/DE_Results_", comp_name, ".csv"),
            row.names = FALSE)
}

# --- 2D. Combined Volcano Figure ---

# Create a more sophisticated layout
figure2 <- wrap_plots(volcano_plots, ncol = 2, nrow = 2) +
  plot_annotation(
    title = 'Differential Expression Analysis',
    subtitle = 'Comparing biTET2/SRSF2 mutations in MDS and CMML (adjusted for age and sex)',
    caption = paste('FDR < ', FDR_THRESHOLD, '; |log2FC| > ', LOGFC_THRESHOLD, sep = ''),
    tag_levels = 'A',
    theme = theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 10)),
      plot.caption = element_text(size = 10, hjust = 1, face = "italic")
    )
  )

ggsave("Publication_Figures/Figure2_Volcano_Plots.pdf",
       figure2, width = 14, height = 12, device = cairo_pdf)
ggsave("Publication_Figures/Figure2_Volcano_Plots.png",
       figure2, width = 14, height = 12, dpi = 300, bg = "white")

cat("\nFigure 2 saved successfully\n")

# --- 2E. Create MA plots as supplementary figures ---

ma_plots <- list()

for (comp_name in names(de_results_list)) {
  res_df <- de_results_list[[comp_name]]
  
  # Create MA plot
  ma_plot <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
    geom_point(aes(color = Significance), size = 1, alpha = 0.6) +
    scale_color_manual(values = sig_colors, name = NULL) +
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    geom_hline(yintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
               linetype = "dashed", color = "grey40", size = 0.5) +
    labs(title = gsub("_", " ", comp_name),
         x = expression(Log[10] ~ "(Mean Expression + 1)"),
         y = expression(Log[2] ~ "Fold Change")) +
    theme_publication() +
    theme(legend.position = "bottom")
  
  ma_plots[[comp_name]] <- ma_plot
}

# Combine MA plots
figure2_supp <- wrap_plots(ma_plots, ncol = 2) +
  plot_annotation(
    title = 'MA Plots - Mean Expression vs Fold Change',
    tag_levels = 'A'
  )

ggsave("Publication_Figures/Figure2_Supp_MA_Plots.pdf",
       figure2_supp, width = 14, height = 10, device = cairo_pdf)

# --- 2F. Summary Statistics ---

# Create summary table of DE results
de_summary <- data.frame(
  Comparison = names(de_results_list),
  Total_Genes_Tested = sapply(de_results_list, nrow),
  Upregulated = sapply(de_results_list, function(x) sum(x$Significance == "Upregulated")),
  Downregulated = sapply(de_results_list, function(x) sum(x$Significance == "Downregulated")),
  Total_DE = sapply(de_results_list, function(x)
    sum(x$Significance %in% c("Upregulated", "Downregulated")))
)

write.csv(de_summary, "Publication_Figures/DE_Summary_Statistics.csv", row.names = FALSE)

cat("\nDifferential Expression Summary:\n")
print(de_summary)
cat("\nAll DE results saved to Publication_Figures/\n")

# --- 2G. Export key objects for downstream analysis ---

# Save the DESeq2 object and results for Part 3 and 4
save(dds, de_results_list, normalized_counts, metadata_de,
     file = "Publication_Figures/DESeq2_results.RData")

cat("\nDESeq2 results saved to DESeq2_results.RData for downstream analysis\n")



################################################################################
### PART 3: PATHWAY ENRICHMENT ANALYSIS (FIXED)
################################################################################

cat("\n--- PART 3: Pathway Enrichment Analysis ---\n")

# Load the DESeq2 results if not already in memory
if (!exists("de_results_list")) {
  load("Publication_Figures/DESeq2_results.RData")
}

# Process each comparison
enrichment_plots <- list()
all_go_results <- list()

for (comp_name in names(de_results_list)) {
  cat("\nProcessing enrichment for:", comp_name, "\n")
  
  results <- de_results_list[[comp_name]]
  
  #rownames(results) = gsub("\\.\\d+$","",rownames(results))
  # Check column names
  cat("  Available columns:", paste(colnames(results), collapse = ", "), "\n")
  
  # Get significant genes - FIX: Use correct column names
  sig_genes <- results %>%
    filter(padj < FDR_THRESHOLD, abs(log2FoldChange) > LOGFC_THRESHOLD) %>%
    pull(Gene)
  
  cat("  Found", length(sig_genes), "significant genes\n")
  
  if (length(sig_genes) < 10) {
    cat("  Too few significant genes for enrichment analysis\n")
    next
  }
  
  # Check if we have Ensembl IDs or gene symbols
  is_ensembl <- any(grepl("^ENSG", sig_genes[1:min(10, length(sig_genes))]))
  
  if (is_ensembl) {
    cat("  Detected Ensembl IDs, converting to gene symbols first\n")
    
    # Clean Ensembl IDs (remove version numbers)
    sig_genes_clean <- gsub("\\.\\d+$", "", sig_genes)
    
    # Convert Ensembl to Symbol first
    symbol_conversion <- tryCatch({
      bitr(sig_genes_clean,
           fromType = "ENSEMBL",
           toType = "SYMBOL",
           OrgDb = org.Hs.eg.db)
    }, error = function(e) {
      cat("  Error in Ensembl to Symbol conversion:", e$message, "\n")
      NULL
    })
    
    if (is.null(symbol_conversion) || nrow(symbol_conversion) == 0) {
      cat("  Failed to convert Ensembl IDs to symbols\n")
      next
    }
    
    # Now convert symbols to Entrez
    entrez_conversion <- tryCatch({
      bitr(symbol_conversion$SYMBOL,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)
    }, error = function(e) {
      cat("  Error in Symbol to Entrez conversion:", e$message, "\n")
      NULL
    })
    
  } else {
    cat("  Detected gene symbols, converting to Entrez IDs\n")
    
    # Direct conversion from Symbol to Entrez
    entrez_conversion <- tryCatch({
      bitr(sig_genes,
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)
    }, error = function(e) {
      cat("  Error in gene ID conversion:", e$message, "\n")
      NULL
    })
  }
  
  if (is.null(entrez_conversion) || nrow(entrez_conversion) == 0) {
    cat("  No genes could be converted to Entrez IDs\n")
    next
  }
  
  cat("  Successfully converted", nrow(entrez_conversion), "genes to Entrez IDs\n")
  
  # GO Enrichment
  cat("  Running GO enrichment analysis...\n")
  go_result <- tryCatch({
    enrichGO(
      gene = entrez_conversion$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
  }, error = function(e) {
    cat("  Error in GO enrichment:", e$message, "\n")
    NULL
  })
  
  if (!is.null(go_result) && nrow(as.data.frame(go_result)) > 0) {
    cat("  Found", nrow(as.data.frame(go_result)), "enriched GO terms\n")
    
    # Store results
    all_go_results[[comp_name]] <- go_result
    
    # Create dot plot
    go_plot <- dotplot(go_result, showCategory = 15) +
      scale_color_viridis_c(name = "Adjusted\nP-value") +
      scale_size_continuous(name = "Gene\nCount", range = c(3, 8)) +
      labs(title = "GO: Biological Process",
           subtitle = gsub("_", " ", comp_name)) +
      theme_publication() +
      theme(axis.text.y = element_text(size = 10))
    
    enrichment_plots[[paste0(comp_name, "_GO")]] <- go_plot
    
    # Save plot
    ggsave(paste0("Publication_Figures/GO_", comp_name, ".pdf"),
           go_plot, width = 10, height = 8, device = cairo_pdf)
    ggsave(paste0("Publication_Figures/GO_", comp_name, ".png"),
           go_plot, width = 10, height = 8, dpi = 300, bg = "white")
    
    # Save results
    go_df <- as.data.frame(go_result)
    write.csv(go_df,
              paste0("Publication_Figures/GO_Results_", comp_name, ".csv"),
              row.names = FALSE)
    
    # Also create a bar plot as alternative visualization
    top_terms <- go_df %>%
      arrange(p.adjust) %>%
      head(20)
    
    bar_plot <- ggplot(top_terms, aes(x = Count, y = reorder(Description, Count))) +
      geom_bar(stat = "identity", fill = cb_palette[1], alpha = 0.8) +
      geom_text(aes(label = Count), hjust = -0.3, size = 3) +
      labs(x = "Gene Count",
           y = "",
           title = "Top Enriched GO Terms",
           subtitle = gsub("_", " ", comp_name)) +
      theme_publication() +
      theme(axis.text.y = element_text(size = 10)) +
      xlim(0, max(top_terms$Count) * 1.1)
    
    ggsave(paste0("Publication_Figures/GO_Bar_", comp_name, ".pdf"),
           bar_plot, width = 10, height = 8, device = cairo_pdf)
    
  } else {
    cat("  No significantly enriched GO terms found\n")
  }
  
  # Try KEGG pathway analysis as well
  cat("  Running KEGG pathway analysis...\n")
  kegg_result <- tryCatch({
    enrichKEGG(
      gene = entrez_conversion$ENTREZID,
      organism = 'hsa',
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    cat("  Error in KEGG enrichment:", e$message, "\n")
    NULL
  })
  
  if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
    cat("  Found", nrow(as.data.frame(kegg_result)), "enriched KEGG pathways\n")
    
    # Create KEGG plot
    kegg_plot <- dotplot(kegg_result, showCategory = 15) +
      scale_color_viridis_c(name = "Adjusted\nP-value") +
      scale_size_continuous(name = "Gene\nCount", range = c(3, 8)) +
      labs(title = "KEGG Pathways",
           subtitle = gsub("_", " ", comp_name)) +
      theme_publication() +
      theme(axis.text.y = element_text(size = 10))
    
    enrichment_plots[[paste0(comp_name, "_KEGG")]] <- kegg_plot
    
    # Save plot
    ggsave(paste0("Publication_Figures/KEGG_", comp_name, ".pdf"),
           kegg_plot, width = 10, height = 8, device = cairo_pdf)
    
    # Save results
    write.csv(as.data.frame(kegg_result),
              paste0("Publication_Figures/KEGG_Results_", comp_name, ".csv"),
              row.names = FALSE)
  }
}

# --- Create combined enrichment figure ---

if (length(enrichment_plots) > 0) {
  # Select GO plots for main figure
  go_plots <- enrichment_plots[grep("_GO$", names(enrichment_plots))]
  
  if (length(go_plots) > 0) {
    figure3 <- wrap_plots(go_plots, ncol = 2) +
      plot_annotation(
        title = 'Gene Ontology Enrichment Analysis',
        subtitle = 'Biological Process terms enriched in differentially expressed genes',
        tag_levels = 'A',
        theme = theme(
          plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
          plot.subtitle = element_text(size = 14, hjust = 0.5)
        )
      )
    
    ggsave("Publication_Figures/Figure3_GO_Enrichment.pdf",
           figure3, width = 16, height = 12, device = cairo_pdf)
    ggsave("Publication_Figures/Figure3_GO_Enrichment.png",
           figure3, width = 16, height = 12, dpi = 300, bg = "white")
    
    cat("\nFigure 3 saved successfully\n")
  }
  
  # Create supplementary figure with KEGG if available
  kegg_plots <- enrichment_plots[grep("_KEGG$", names(enrichment_plots))]
  
  if (length(kegg_plots) > 0) {
    figure3_supp <- wrap_plots(kegg_plots, ncol = 2) +
      plot_annotation(
        title = 'KEGG Pathway Enrichment Analysis',
        tag_levels = 'A'
      )
    
    ggsave("Publication_Figures/Figure3_Supp_KEGG_Enrichment.pdf",
           figure3_supp, width = 16, height = 12, device = cairo_pdf)
  }
}

# --- Create summary of enrichment results ---

enrichment_summary <- data.frame(
  Comparison = character(),
  DE_Genes = numeric(),
  Converted_to_Entrez = numeric(),
  GO_Terms_Enriched = numeric(),
  KEGG_Pathways_Enriched = numeric()
)

for (comp_name in names(de_results_list)) {
  results <- de_results_list[[comp_name]]
  
  # Count DE genes
  n_de <- sum(results$padj < FDR_THRESHOLD &
              abs(results$log2FoldChange) > LOGFC_THRESHOLD, na.rm = TRUE)
  
  # Get GO results
  go_count <- 0
  if (comp_name %in% names(all_go_results)) {
    go_count <- nrow(as.data.frame(all_go_results[[comp_name]]))
  }
  
  # Get KEGG count (if saved)
  kegg_count <- 0
  kegg_file <- paste0("Publication_Figures/KEGG_Results_", comp_name, ".csv")
  if (file.exists(kegg_file)) {
    kegg_data <- read.csv(kegg_file)
    kegg_count <- nrow(kegg_data)
  }
  
  # Get conversion count
  conv_count <- 0
  go_file <- paste0("Publication_Figures/GO_Results_", comp_name, ".csv")
  if (file.exists(go_file)) {
    go_data <- read.csv(go_file)
    if (nrow(go_data) > 0) {
      # Estimate from gene ratio
      conv_count <- round(mean(as.numeric(sub("/.*", "", go_data$GeneRatio))))
    }
  }
  
  enrichment_summary <- rbind(enrichment_summary, data.frame(
    Comparison = comp_name,
    DE_Genes = n_de,
    Converted_to_Entrez = conv_count,
    GO_Terms_Enriched = go_count,
    KEGG_Pathways_Enriched = kegg_count
  ))
}

write.csv(enrichment_summary,
          "Publication_Figures/Enrichment_Summary.csv",
          row.names = FALSE)

cat("\nEnrichment Analysis Summary:\n")
print(enrichment_summary)
cat("\nAll enrichment results saved to Publication_Figures/\n")

cat("\nAll analyses completed successfully!\n")
