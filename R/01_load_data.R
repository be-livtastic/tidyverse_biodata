# ============================================================
# 01_load_data.R
# Purpose: Load and inspect the gene expression dataset
# ============================================================

library(tidyverse)

# ----- Load data -----
data_path <- here::here("data", "gene_expression.csv")

gene_expr <- read_csv(data_path, show_col_types = FALSE)

# ----- Initial inspection -----
glimpse(gene_expr)
cat("\nDimensions:", nrow(gene_expr), "genes x", ncol(gene_expr), "columns\n")

# Column names
cat("\nColumn names:\n")
print(colnames(gene_expr))

# First few rows
cat("\nFirst 6 rows:\n")
print(head(gene_expr))

# Check for missing values
cat("\nMissing values per column:\n")
print(colSums(is.na(gene_expr)))

# Gene biotype summary
cat("\nGene biotype counts:\n")
print(count(gene_expr, gene_biotype))

# ----- Save as RDS for downstream scripts -----
saveRDS(gene_expr, here::here("data", "gene_expression.rds"))
cat("\nData saved to data/gene_expression.rds\n")
