install.packages("tidyverse")
install.packages("pheatmap")

library(tidyverse)
library(pheatmap)
library(ggplot2)
library(dplyr)

### Data Wrangling

# Load expression data
expr_data <- read.csv("week2_messy_expression_data.csv")

# Load metadata
metadata <- read.csv("week2_sample_metadata.csv")

glimpse(expr_data)
glimpse(metadata)

### Dataset Overview
# Number of genes (rows) and samples (columns)
cat("Genes:", nrow(expr_data), "\n")
cat("Samples:", ncol(expr_data) - 1, "\n") # subtract 1 if first col is gene names # nolint: line_length_linter.

# Samples (columns) with any NAs
na_per_sample <- colSums(is.na(expr_data))
cat("Samples with NAs:", sum(na_per_sample > 0), "\n")
print(na_per_sample[na_per_sample > 0])

# Range of expression values (excluding non-numeric gene name column)
expr_matrix <- expr_data[, sapply(expr_data, is.numeric)]
cat("Expression range:", range(expr_matrix, na.rm = TRUE), "\n")


