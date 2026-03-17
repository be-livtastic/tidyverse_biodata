# ============================================================
# 02_wrangle_data.R
# Purpose: Reshape and enrich the gene expression dataset
#          using tidyverse tools
# ============================================================

library(tidyverse)

# ----- Load data -----
gene_expr <- readRDS(here::here("data", "gene_expression.rds"))

# ----- Pivot to long (tidy) format -----
gene_expr_long <- gene_expr |>
  pivot_longer(
    cols = ctrl_rep1:treat_rep3,
    names_to  = "sample",
    values_to = "expression"
  ) |>
  # Extract condition and replicate from the sample column
  mutate(
    condition  = if_else(str_starts(sample, "ctrl"), "Control", "Treatment"),
    replicate  = str_extract(sample, "\\d+$"),
    .after     = sample
  )

cat("Long-format dimensions:", nrow(gene_expr_long), "rows x",
    ncol(gene_expr_long), "columns\n")
print(head(gene_expr_long))

# ----- Compute per-gene mean expression per condition -----
gene_means <- gene_expr_long |>
  group_by(gene_id, gene_name, gene_biotype, condition) |>
  summarise(mean_expr = mean(expression), .groups = "drop")

# ----- Compute log2 fold change (Treatment vs Control) -----
gene_fc <- gene_means |>
  pivot_wider(names_from = condition, values_from = mean_expr) |>
  mutate(
    log2fc       = Treatment - Control,          # already log2-scale data
    abs_log2fc   = abs(log2fc),
    direction    = case_when(
      log2fc >  1 ~ "Up",
      log2fc < -1 ~ "Down",
      TRUE        ~ "Stable"
    )
  )

cat("\nFold-change table (first 6 rows):\n")
print(head(gene_fc))

cat("\nDEG direction counts:\n")
print(count(gene_fc, direction))

# ----- Save wrangled objects -----
saveRDS(gene_expr_long, here::here("data", "gene_expression_long.rds"))
saveRDS(gene_fc,        here::here("data", "gene_fc.rds"))
cat("\nWrangled data saved.\n")
