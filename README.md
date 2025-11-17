# Belief Updating Project

> Computational modeling of human belief updating using Approximate Bayesian Computation and Bayesian Model Selection

[![R](https://img.shields.io/badge/R-%E2%89%A54.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## Overview

This project compares **17 cognitive models** of belief updating across:
- Two datasets: experimental and replication study
- Two presentation formats: probability and frequency

### Methodology

- **ABC-RF**: Random Forest for individual-level model inference
- **ABC**: Rejection sampling for parameter estimation  
- **RFX-BMS**: Random Effects Bayesian Model Selection

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

### 2. Run Complete Pipeline

```r
setwd("/path/to/belief-update/")

# Complete analysis workflow
source("scripts/00_data_process.R")                 # Clean data
source("scripts/01_simulate_parameters.R")          # Generate priors
source("scripts/02_generate_predictions.R")         # Model predictions
source("scripts/03_summary_stats.R")                # Summary statistics
source("scripts/04_abcrf_individual_inference.R")   # ABC-RF Individual Model Fitting
source("scripts/05_posterior.R")                    # Posterior Check
source("scripts/05_posterior_MH.R")                 # Posterior Check

# Generate outputs
source("scripts/00_plot.R")                         # Generate all figures
source("scripts/00_stats_table.R")                  # Generate all tables
```

### 3. Quick Results (Pre-computed)

If `.rds` files exist in `data/`, skip steps 1-7 and directly generate outputs:

```r
source("scripts/00_plot.R")         # Generate figures
source("scripts/00_stats_table.R")  # Generate tables
```

---

## Cognitive Models

### Summary of Candidate Models

<details>
<summary><b>Click to expand full model table</b></summary>

| **Family** | **Model / Category** | **Abbrev.** | **Description** |
|------------|---------------------|-------------|-----------------|
| **Bayes Rule** | Bayes Rule | `Bayes` | The normative Bayesian posterior, serving as the optimal benchmark for comparison across models. |
| **Simple Heuristics** | Heuristic rules | `BO`, `REP`, `FC` | Simple rules: Base rate Only (BO), REPresent (REP), False alarm Complement (FC). |
| **Complex Heuristics** | Lexicographic | `LE` | Sequential cue evaluation with thresholds on BR, HR, and FAR. |
| | Adaptive heuristic | `AH` | Selects, on each trial, the heuristic that minimizes absolute deviation from the Bayesian posterior. |
| | Mixture heuristic | `MH` | Draws a heuristic on each trial from a categorical mixture with Dirichlet-weighted probabilities. |
| **Linear Additive** | Linear additive | `LA` | Weighted linear combination of BR, HR, and FAR cues. |
| **Bayesian Sampler** | Symmetric prior | `BS` | Limited-sample Bayesian inference with a symmetric prior centered at 0.5. |
| | Asymmetric prior | `BS-A` | Limited-sample Bayesian inference with a prior mean drawn from a flexible subjective distribution. |
| **Heuristic-Anchored BS** | Single-anchor | `HABS_anchor` | Dual-process mixture: the response is either the heuristic anchor or an anchored Bayesian sample. The subscript *anchor* denotes a specific heuristic (e.g., BO, REP, FC, JO, LS, or 0.5). |
| | Mixture-anchor | `HABS_mixed` | Extension of HABS_anchor in which the anchor α is drawn on each trial from a set of simple heuristics according to probability vector **p**. Captures individual differences in deliberation (q) and anchor preferences (**p**). |

</details>

### Simple Heuristics Reference

| Abbrev. | Full Name | Formula/Description |
|---------|-----------|---------------------|
| `HO` | Hit rate Only | Response = HR |
| `FO` | False alarm complement | Response = 1 - FAR |
| `JO` | Joint Occurrence | Response = BR × HR |
| `LS` | Likelihood Subtraction | Response = HR - FAR |
| `H` | 50% Heuristic | Always responds 0.5 |

---

## Analysis Scripts

- `00_data_process.R`
- `01_simulate_parameters.R`
- `02_generate_predictions.R`
- `03_summary_stats.R`
- `04_abcrf_individual_inference.R`
- `05_posterior.R`
- `05_posterior_MH.R`
- `00_plot.R`
- `00_stats_table.R`

---

## Outputs

| Type | Location | Format |
|------|----------|--------|
| **Figures (high-res)** | `fig/tiff/` | TIFF (600 DPI) |
| **Figures (web)** | `fig/png/` | PNG (300 DPI) |
| **Tables** | `Table/*.tex` | LaTeX format |
| **Data** | `data/*.rds` | R data files |

---

## System Requirements

| Component | Requirement |
|-----------|-------------|
| **R Version** | ≥ 4.0 |
| **RAM** | 16 GB minimum (32 GB recommended) |
| **CPU** | Multi-core processor (uses parallelization) |
| **Disk Space** | ~5 GB |

---

## Author

**Yitong Lin**  
Contact: yitong.lin@warwick.ac.uk

---

## Troubleshooting

<details>
<summary><b>Common Issues and Solutions</b></summary>

### Out of Memory Errors

**Problem**: R crashes with memory allocation errors  
**Solution**: Reduce the number of cores used in parallelized scripts

```r
# In scripts, modify:
n_cores <- detectCores() - 2  # Use fewer cores
```

### Slow Performance

**Problem**: Analysis pipeline takes excessively long  
**Solution**: Use pre-computed `.rds` files from the `data/` directory

### Missing Package Dependencies

**Problem**: Package installation fails  
**Solution**: Install packages with full dependencies

```r
install.packages("package_name", dependencies = TRUE)
```

### Large File Warnings (Git)

**Problem**: Git warns about large TIFF files (>50 MB)  
**Solution**: Consider using Git Large File Storage (LFS)

```bash
git lfs install
git lfs track "*.tiff"
git add .gitattributes
```

</details>

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

---

<div align="center">

**Belief Updating Project**  
Computational Cognitive Science Laboratory

</div>
