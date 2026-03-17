# ============================================================
# 02_wrangle_data.R
# Purpose: Reshape and enrich the gene expression dataset
#          using tidyverse tools
# ============================================================

library(tidyverse)

### Filter Low Quality Samples
high_quality_samples <- metadata$sample_id[metadata$quality_score >= 0.85]
high_quality_samples

# Validate and align metadata sample IDs with expression columns.
expr_sample_cols <- names(expr_data)[-1]
expr_sample_cols
expr_sample_cols_norm <- trimws(as.character(expr_sample_cols))
high_quality_samples_norm <- trimws(as.character(high_quality_samples))

matched_samples <- high_quality_samples[high_quality_samples_norm %in% expr_sample_cols_norm]
missing_samples <- high_quality_samples[!(high_quality_samples_norm %in% expr_sample_cols_norm)]

if (length(missing_samples) > 0) {
    warning(
        paste0(
            "These high-quality metadata sample IDs are missing in expr_data columns: ", # nolint
            paste(unique(missing_samples), collapse = ", ")
        ) # nolint: indentation_linter.
    )
}

if (length(matched_samples) == 0) {
    stop(
        "No high-quality metadata sample IDs matched expr_data column names. Check metadata$sample_id formatting and expression column names." # nolint
    )
}

# Keep gene column + only matched high-quality sample columns
expr_data_filtered <- expr_data[, c(names(expr_data)[1], matched_samples), drop = FALSE]
expr_data_filtered

cat("High-quality samples requested:", length(high_quality_samples), "\n")
cat("High-quality samples matched:", length(matched_samples), "\n")
cat("High-quality samples missing:", length(unique(missing_samples)), "\n")

cat("Samples before filtering:", ncol(expr_data) - 1, "\n")
cat("Samples after filtering:", ncol(expr_data_filtered) - 1, "\n")

### Per-gene Summary Statistics
gene_col <- names(expr_data_filtered)[1] # Assuming the first column is gene names
gene_col

# Pivot to long format and calculate mean, variance, and NA count per gene
gene_pivot <- expr_data_filtered %>%
    pivot_longer(cols = -all_of(gene_col), names_to = "sample", values_to = "expression") %>% #exclude gene name column from pivot
    group_by(across(all_of(gene_col)))
gene_pivot

gene_stats <- gene_pivot %>%
    summarize(
        mean_expr = mean(expression, na.rm = TRUE),
        variance  = var(expression, na.rm = TRUE),
        n_na      = sum(is.na(expression)), #number of NAs per gene
        .groups   = "drop" #to ungroup after summarization
    )

print(gene_stats)
