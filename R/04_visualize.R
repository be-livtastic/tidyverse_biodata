# ============================================================
# 04_visualize.R
# Purpose: Generate ggplot2 visualisations for gene expression
#          data, including boxplots, a volcano plot, a heatmap,
#          and a PCA plot
# ============================================================

library(tidyverse)

# Create output directory if it does not already exist
if (!dir.exists(here::here("output"))) dir.create(here::here("output"))

# ----- Load data -----
gene_expr_long <- readRDS(here::here("data", "gene_expression_long.rds"))
gene_fc        <- readRDS(here::here("data", "gene_fc.rds"))
gene_expr_wide <- readRDS(here::here("data", "gene_expression.rds"))

# ── 1. Boxplot: expression distribution per sample ────────────────────────────
p_box <- ggplot(gene_expr_long,
                aes(x = sample, y = expression, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = c(Control = "#4393C3", Treatment = "#D6604D")) +
  labs(
    title    = "Gene Expression Distribution per Sample",
    subtitle = "Log2-normalised expression across all 50 genes",
    x        = "Sample",
    y        = "Expression (log2)",
    fill     = "Condition"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here::here("output", "01_boxplot_samples.png"),
       plot = p_box, width = 8, height = 5, dpi = 150)
cat("Saved: 01_boxplot_samples.png\n")

# ── 2. Volcano plot ───────────────────────────────────────────────────────────
volcano_data <- gene_fc |>
  mutate(
    label = if_else(abs_log2fc >= 3, gene_name, NA_character_)
  )

p_volcano <- ggplot(volcano_data,
                    aes(x = log2fc, y = abs_log2fc,
                        colour = direction, label = label)) +
  geom_point(alpha = 0.8, size = 2.5) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 15, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(
    values = c(Up = "#D6604D", Down = "#4393C3", Stable = "grey60")
  ) +
  labs(
    title    = "Volcano Plot: Treatment vs Control",
    x        = "Log2 Fold Change",
    y        = "|Log2 Fold Change|",
    colour   = "Direction"
  ) +
  theme_bw(base_size = 12)

ggsave(here::here("output", "02_volcano_plot.png"),
       plot = p_volcano, width = 7, height = 5, dpi = 150)
cat("Saved: 02_volcano_plot.png\n")

# ── 3. Heatmap of top 20 variable genes ──────────────────────────────────────
top20_genes <- gene_expr_long |>
  group_by(gene_id) |>
  summarise(cv = sd(expression) / mean(expression), .groups = "drop") |>
  slice_max(cv, n = 20) |>
  pull(gene_id)

heatmap_data <- gene_expr_wide |>
  filter(gene_id %in% top20_genes) |>
  pivot_longer(ctrl_rep1:treat_rep3,
               names_to = "sample", values_to = "expression")

p_heat <- ggplot(heatmap_data,
                 aes(x = sample, y = fct_reorder(gene_name, expression),
                     fill = expression)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "#4393C3",
    mid      = "white",
    high     = "#D6604D",
    midpoint = median(heatmap_data$expression)
  ) +
  labs(
    title = "Heatmap: Top 20 Most Variable Genes",
    x     = "Sample",
    y     = "Gene",
    fill  = "Expression\n(log2)"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here::here("output", "03_heatmap_top20.png"),
       plot = p_heat, width = 7, height = 7, dpi = 150)
cat("Saved: 03_heatmap_top20.png\n")

# ── 4. PCA plot ───────────────────────────────────────────────────────────────
# Build a samples × genes matrix for PCA
pca_matrix <- gene_expr_wide |>
  select(gene_id, ctrl_rep1:treat_rep3) |>
  column_to_rownames("gene_id") |>
  t()                              # samples as rows, genes as columns

pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)

pca_df <- as_tibble(pca_result$x[, 1:2], rownames = "sample") |>
  mutate(condition = if_else(str_starts(sample, "ctrl"), "Control", "Treatment"))

# Variance explained
var_exp <- summary(pca_result)$importance["Proportion of Variance", ] * 100

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                             colour = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3.5) +
  scale_colour_manual(values = c(Control = "#4393C3", Treatment = "#D6604D")) +
  labs(
    title  = "PCA of Gene Expression Profiles",
    x      = sprintf("PC1 (%.1f%% variance)", var_exp[1]),
    y      = sprintf("PC2 (%.1f%% variance)", var_exp[2]),
    colour = "Condition"
  ) +
  theme_bw(base_size = 12)

ggsave(here::here("output", "04_pca_plot.png"),
       plot = p_pca, width = 6, height = 5, dpi = 150)
cat("Saved: 04_pca_plot.png\n")

# ── 5. Bar chart: DEG counts by biotype ──────────────────────────────────────
p_bar <- gene_fc |>
  filter(direction != "Stable") |>
  count(direction, gene_biotype) |>
  ggplot(aes(x = gene_biotype, y = n, fill = direction)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c(Up = "#D6604D", Down = "#4393C3")) +
  labs(
    title = "Differentially Expressed Genes by Biotype",
    x     = "Gene Biotype",
    y     = "Number of DEGs",
    fill  = "Direction"
  ) +
  theme_bw(base_size = 12)

ggsave(here::here("output", "05_deg_bar_biotype.png"),
       plot = p_bar, width = 6, height = 4, dpi = 150)
cat("Saved: 05_deg_bar_biotype.png\n")

cat("\nAll plots saved to output/\n")
