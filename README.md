# tidyverse_biodata

A hands-on project for **loading, wrangling, summarising, and visualising gene expression data** using the [tidyverse](https://www.tidyverse.org/) ecosystem in R.

---

## Project structure

```
tidyverse_biodata/
├── tidyverse_biodata.Rproj   # RStudio project file
├── data/
│   └── gene_expression.csv   # Simulated log2-normalised RNA-seq dataset
├── R/
│   ├── 01_load_data.R        # Load CSV, inspect, save as RDS
│   ├── 02_wrangle_data.R     # Pivot to long format, compute fold changes
│   ├── 03_summary_stats.R    # Per-sample and per-gene summary statistics
│   └── 04_visualize.R        # ggplot2 figures (boxplot, volcano, heatmap, PCA)
├── reports/
│   └── gene_expression_report.Rmd  # Self-contained R Markdown report
└── output/                   # Generated CSV tables and PNG figures (git-ignored)
```

---

## Dataset

`data/gene_expression.csv` contains **50 genes** (47 protein-coding, 3 lncRNA)
across **6 samples**: three replicates of a *Control* condition
(`ctrl_rep1–3`) and three replicates of a *Treatment* condition
(`treat_rep1–3`). Expression values are already on the **log₂ scale**.

| Column | Description |
|---|---|
| `gene_id` | Ensembl-style gene identifier |
| `gene_name` | HGNC gene symbol |
| `gene_biotype` | `protein_coding` or `lncRNA` |
| `ctrl_rep1–3` | Log2 expression in control replicates |
| `treat_rep1–3` | Log2 expression in treatment replicates |

---

## Scripts

Run the scripts in order from the project root (or open the `.Rproj` file and
use the RStudio console/Build panel):

```r
source("R/01_load_data.R")     # saves data/gene_expression.rds
source("R/02_wrangle_data.R")  # saves data/gene_expression_long.rds & gene_fc.rds
source("R/03_summary_stats.R") # writes CSV tables to output/
source("R/04_visualize.R")     # writes PNG figures to output/
```

### What each script does

| Script | Key tidyverse functions used |
|---|---|
| `01_load_data.R` | `read_csv()`, `glimpse()`, `count()` |
| `02_wrangle_data.R` | `pivot_longer()`, `mutate()`, `str_starts()`, `group_by()`, `summarise()`, `pivot_wider()`, `case_when()` |
| `03_summary_stats.R` | `group_by()`, `summarise()`, `slice_max()`, `slice_min()`, `write_csv()` |
| `04_visualize.R` | `ggplot2`, `geom_boxplot()`, `geom_tile()`, `geom_point()`, `geom_col()`, `scale_fill_gradient2()`, `prcomp()`, `ggsave()` |

---

## R Markdown report

Open `reports/gene_expression_report.Rmd` in RStudio and click **Knit** to
produce a self-contained HTML report that walks through every analysis step
with narrative, code, tables, and figures.

---

## Requirements

```r
install.packages(c("tidyverse", "ggrepel", "here", "knitr", "rmarkdown"))
```

- R ≥ 4.1  
- tidyverse ≥ 2.0  
- ggrepel ≥ 0.9  
- here ≥ 1.0  

---

## Output files

After running all scripts the `output/` directory will contain:

| File | Description |
|---|---|
| `sample_summary.csv` | Mean, median, SD per sample |
| `gene_summary.csv` | CV-ranked gene variability table |
| `gene_foldchange.csv` | Log2FC, direction for all genes |
| `01_boxplot_samples.png` | Expression distribution boxplots |
| `02_volcano_plot.png` | Volcano plot (Treatment vs Control) |
| `03_heatmap_top20.png` | Heatmap of top 20 variable genes |
| `04_pca_plot.png` | PCA of sample expression profiles |
| `05_deg_bar_biotype.png` | DEG counts by gene biotype |

---

## License

[MIT](LICENSE)
