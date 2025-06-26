# Belief Updating Project

This repository contains R code and data for simulating, analyzing, and visualizing cognitive models of belief updating.

## Structure

- `data/`: raw and processed datasets (e.g., df.xlsx)
- `R/`: model definitions and evaluation metrics
- `scripts/`: main analysis pipeline
- `figs/`: generated figures
- `tests/`: unit tests with testthat
- `.gitignore`: standard ignore rules
- `belief-update.Rproj`: RStudio project file

## Setup

\```r
# Install dependencies
install.packages(c('tidyverse', 'data.table', 'ggplot2', 'testthat'))

# (Optional) Initialize renv for reproducibility
renv::init()
\```

## Usage

Run each script in `scripts/` in order:

1. `01_simulate_parameters.R`
2. `02_generate_predictions.R`
3. `03_analysis_plots.R`

## License

MIT License. Attribution appreciated.
