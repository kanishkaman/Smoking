# Gene Expression Analysis Using 2-Way ANOVA

This repository contains code and analysis for exploring the interaction effects of smoking status and gender on gene expression. Using a two-way ANOVA, we identify genes that may respond differently to smoking in males and females.

## Repository Contents

- **`assignment3.py`**: Main script for loading data, performing two-way ANOVA, and visualizing p-values.
  - Loads and preprocesses data (reverts log transformation).
  - Defines design matrices for additive and interaction models.
  - Calculates F-statistics and p-values for each gene.
- **`data/`**: Contains the data file.
  - `Raw Data_GeneSpring.txt`: The gene expression dataset used in the analysis.
- **`Plots/`**: Contains histogram plots of interaction p-values.
  - `Plot1.png`: Using preprocessed (exponentiated) values.
  - `Plot2.png`: Using original log-transformed values.
- **`Report3.pdf`**: Full report on methodology, data, and results.

## Quick Start

1. Ensure the `Raw Data_GeneSpring.txt` file is located in the `data/` folder.
2. Run `assignment3.py` to execute the analysis and generate histograms.

This project highlights how smoking and gender may influence gene expression through interaction effects, providing insights into biological responses.

