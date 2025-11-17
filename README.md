# Belief Updating Project

Computational modeling of human belief updating using Approximate Bayesian Computation and Bayesian Model Selection.

## Overview

This project compares 17 cognitive models of belief updating across two datasets (experimental and replication) and two presentation formats (probability and frequency).

**Methods:**
- ABC Random Forest for individual-level model inference
- ABC with rejection sampling for parameter estimation
- Random Effects Bayesian Model Selection

---

## Project Structure

```
belief-update/
├── data/              # Raw and processed data (.csv, .rds)
├── R/                 # Model definitions and helper functions
├── scripts/           # Analysis pipeline (run in order)
├── fig/               # Generated figures (TIFF, PNG)
└── Table/             # LaTeX tables
```

---

## Quick Start

### 1. Install Dependencies

```r
install.packages("pacman")
pacman::p_load(
  data.table, tidyverse, abc, abcrf, ranger, MCMCpack,
  ggplot2, patchwork, parallel, doParallel, knitr, kableExtra
)
```

### 2. Run Analysis Pipeline

```r
setwd("/path/to/belief-update/")

# Complete workflow (8-10 hours)
source("scripts/00_data_process.R")         # Clean data
source("scripts/01_simulate_parameters.R")  # Generate priors
source("scripts/02_generate_predictions.R") # Model predictions
source("scripts/03_summary_stats.R")        # Summary statistics
source("scripts/04_abcrf_individual_inference.R")  # ABC-RF
source("scripts/05_posterior.R")            # ABC inference (HAmix)
source("scripts/05_posterior_MH.R")         # ABC-MH (MH model)

# Generate outputs
source("scripts/00_plot.R")                 # All figures
source("scripts/00_stats_table.R")          # All tables
```

### 3. Quick Start (Pre-computed Results)

If `.rds` files exist in `data/`:

```r
source("scripts/00_plot.R")         # Generate figures
source("scripts/00_stats_table.R")  # Generate tables
```

---

## Models

### Summary of Candidate Models

Models are grouped by family with abbreviations and formal descriptions.

| **Family** | **Model / Category** | **Abbrev.** | **Description** |
|------------|---------------------|-------------|-----------------|
| **Bayes Rule** | Bayes Rule | Bayes | The normative Bayesian posterior, serving as the optimal benchmark for comparison across models. |
| **Simple Heuristics** | Heuristic rules | BO, REP, FC | Simple rules: Base rate Only (BO), REPresent (REP), False alarm Complement (FC). |
| **Complex Heuristics** | Lexicographic | LE | Sequential cue evaluation with thresholds on BR, HR, and FAR. |
| | Adaptive heuristic | AH | Selects, on each trial, the heuristic that minimizes absolute deviation from the Bayesian posterior. |
| | Mixture heuristic | MH | Draws a heuristic on each trial from a categorical mixture with Dirichlet-weighted probabilities. |
| **Linear Additive** | Linear additive | LA | Weighted linear combination of BR, HR, and FAR cues. |
| **Bayesian Sampler** | Symmetric prior | BS | Limited-sample Bayesian inference with a symmetric prior centered at 0.5. |
| | Asymmetric prior | BS-A | Limited-sample Bayesian inference with a prior mean drawn from a flexible subjective distribution. |
| **Heuristic-Anchored Bayesian Sampler (HABS)** | Single-anchor | HABS_anchor | Dual-process mixture: the response is either the heuristic anchor or an anchored Bayesian sample. The subscript *anchor* denotes a specific heuristic (e.g., BO, REP, FC, JO, LS, or 0.5). |
| | Mixture-anchor | HABS_mixed | Extension of HABS_anchor in which the anchor α is drawn on each trial from a set of simple heuristics according to probability vector **p**. Captures individual differences in deliberation (q) and anchor preferences (**p**). |

### Additional Heuristics

- **HO**: Hit rate Only
- **FO**: False alarm rate complement (1 - FAR)
- **JO**: Joint Occurrence (BR × HR)
- **LS**: Likelihood Subtraction (HR - FAR)
- **H**: 50% Heuristic (always responds 0.5)

---

## Key Scripts

| Script | Purpose | Runtime |
|--------|---------|---------|
| `01_simulate_parameters.R` | Generate 50K prior samples | ~30 min |
| `02_generate_predictions.R` | Compute model predictions | ~2 hours |
| `03_summary_stats.R` | Calculate ABC features | ~10 min |
| `04_abcrf_individual_inference.R` | ABC-RF classification | ~1 hour |
| `05_posterior.R` | ABC parameter estimation | ~2 hours |
| `05_posterior_MH.R` | ABC-MH inference | ~1 hour |
| `00_plot.R` | Generate all figures | ~5 min |
| `00_stats_table.R` | Generate LaTeX tables | ~2 min |

---

## Outputs

- **Figures**: `fig/tiff/` (600 DPI), `fig/png/` (300 DPI)
- **Tables**: `Table/*.tex` (LaTeX format)
- **Data**: `data/*.rds` (posteriors, parameters, predictions)

---

## System Requirements

- **R**: ≥4.0
- **RAM**: 16 GB minimum (32 GB recommended)
- **CPU**: Multi-core (scripts use parallelization)
- **Disk**: ~5 GB

---

## Citation

```
[Your publication details]
```

---

## Author

Yitong Lin | Created: June 26, 2025

---

## Troubleshooting

**Out of memory?** Reduce `n_cores` in scripts  
**Slow performance?** Use pre-computed `.rds` files  
**Missing packages?** Run `pacman::p_load()` with `dependencies = TRUE`

For issues: [your.email@domain.com]
