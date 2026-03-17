
**Tidyverse BioData Wrangling: Learning R with Gene Expression Data**

---

# Description

This project was developed as part of my learning journey in tidyverse-based data wrangling and biological data analysis. The focus is on building clean, reproducible workflows for structured research datasets.

This repository demonstrates how messy biological datasets can be cleaned, reshaped, and summarized to extract meaningful insights. The workflow includes inspecting data quality, transforming wide datasets into tidy format, and computing summary statistics for genes across experimental conditions.

The dataset used here is **AI-generated to simulate a transcriptomics experiment**, allowing safe experimentation with data cleaning and biological interpretation.

The project highlights practical skills in:

* data wrangling
* exploratory data analysis
* tidy data principles
* biological data interpretation
* reproducible research workflows.

---

# Project Goals

This repository was created as a learning exercise to:

* Practice using the **tidyverse workflow in R**
* Understand how **biological datasets are structured**
* Apply **data wrangling techniques to messy experimental data**
* Generate interpretable summary statistics for gene expression
* Build a clear and reproducible project structure suitable for collaborative research.

---

# Technologies Used

This project uses **R** and several packages from the **tidyverse**.

Key packages include:

* **dplyr** – data manipulation and summarization
* **tidyr** – reshaping datasets
* **ggplot2** – visualization
* **pheatmap** – heatmap generation for gene expression data.

---

# Dataset Overview

The repository uses an **AI-generated gene expression dataset** simulating a small transcriptomics experiment.

Expression data includes intentionally messy elements such as:

* missing values
* uneven sample sizes
* variability across experimental conditions.

This allows practice with **realistic data cleaning challenges** commonly encountered in biological research.

---

# Example Results

Example gene summary output:

| gene_id  | mean_expr | variance | n_na |
| -------- | --------- | -------- | ---- |
| GENE_001 | 592       | 165253   | 0    |
| GENE_003 | 1274      | 1181     | 0    |
| GENE_007 | 2.08      | 0.022    | 0    |

Interpretation:

* **GENE_003** shows very high expression across samples
* **GENE_001** has extremely high variance, suggesting strong changes between control and treatment samples
* **GENE_007** shows very low expression and minimal variability.

These patterns mimic typical biological behaviors seen in transcriptomics experiments.

---

# Key Learning Outcomes

Through this project I practiced:

### Data Wrangling

* filtering datasets
* reshaping data
* grouping and summarizing observations.

### Tidy Data Principles

Understanding how to structure data so each:

* variable forms a column
* observation forms a row
* value occupies a single cell.

### Biological Data Interpretation

Learning how summary statistics can reveal:

* differential expression
* gene stability
* potential experimental effects.

### Reproducible Research

Organizing analysis scripts and data files in a clear project structure.


---

# Who This Repository Is For

This project may be useful for:

* students learning **R for bioinformatics**
* beginners exploring **tidyverse workflows**
* researchers interested in **data wrangling biological datasets**
* anyone curious about **reproducible computational biology practices**.
