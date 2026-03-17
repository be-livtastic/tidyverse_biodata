# ============================================================
# 03_summary_stats.R
# Purpose: Compute summary statistics for the gene expression
#          dataset (per-gene and per-sample)
# ============================================================

library(tidyverse)

# ----- Load data -----
gene_expr_long <- readRDS(here::here("data", "gene_expression_long.rds"))
gene_fc        <- readRDS(here::here("data", "gene_fc.rds"))

# ----- Per-sample summary -----
sample_summary <- gene_expr_long |>
  group_by(sample, condition) |>
  summarise(
    n_genes  = n(),
    mean     = mean(expression),
    median   = median(expression),
    sd       = sd(expression),
    min      = min(expression),
    max      = max(expression),
    .groups  = "drop"
  ) |>
  arrange(condition, sample)

cat("=== Per-sample summary statistics ===\n")
print(sample_summary)

# ----- Per-gene summary (across all samples) -----
gene_summary <- gene_expr_long |>
  group_by(gene_id, gene_name, gene_biotype) |>
  summarise(
    mean_all  = mean(expression),
    sd_all    = sd(expression),
    cv        = sd_all / mean_all,   # coefficient of variation
    .groups   = "drop"
  ) |>
  arrange(desc(cv))

cat("\n=== Top 10 most variable genes (by CV) ===\n")
print(head(gene_summary, 10))

# ----- Differentially expressed gene (DEG) summary -----
deg_summary <- gene_fc |>
  group_by(direction, gene_biotype) |>
  summarise(n = n(), .groups = "drop") |>
  arrange(direction, gene_biotype)

cat("\n=== DEG counts by direction and biotype ===\n")
print(deg_summary)

# ----- Top upregulated and downregulated genes -----
top_up <- gene_fc |>
  filter(direction == "Up") |>
  slice_max(log2fc, n = 10) |>
  select(gene_name, gene_biotype, Control, Treatment, log2fc)

top_down <- gene_fc |>
  filter(direction == "Down") |>
  slice_min(log2fc, n = 10) |>
  select(gene_name, gene_biotype, Control, Treatment, log2fc)

cat("\n=== Top 10 upregulated genes ===\n")
print(top_up)

cat("\n=== Top 10 downregulated genes ===\n")
print(top_down)

# ----- Save summaries -----
write_csv(sample_summary, here::here("output", "sample_summary.csv"))
write_csv(gene_summary,   here::here("output", "gene_summary.csv"))
write_csv(gene_fc,        here::here("output", "gene_foldchange.csv"))

cat("\nSummary tables saved to output/\n")
