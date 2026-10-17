# Integrating Heuristics and Mental Sampling in Belief Updating

This repository contains the data analysis and computational modelling code for the research project **Integrating Heuristics and Mental Sampling in Belief Updating**, submitted to *Journal of Experimental Psychology: General*.

## Overview

This project investigates how people combine simple heuristics and Bayesian sampling when updating beliefs in response to new evidence.

We examine whether probability judgements can be explained by a hybrid computational framework in which people use heuristics to generate initial estimates and refine them through limited Bayesian sampling.

The repository includes analyses from three datasets:

1. **Stengård dataset reanalysis:** Reanalysis of an existing dataset to evaluate heuristic, Bayesian, sampling, and hybrid models.
2. **Experiment 1:** An experimental investigation of belief updating using frequency and probability presentation formats.
3. **Experiment 2:** A follow-up experiment examining how time pressure influences probability judgements, heuristic use, and Bayesian sampling.

Supporting descriptive analyses also use `data/sirota.csv` (Sirota et al., 2024), particularly Figure 1. The three principal datasets and supporting input are mapped in the [analysis guide](docs/analysis_pipeline.md#study-mapping-and-glossary).

## Repository Structure

```text
R/                         Shared model, inference, statistics, and reporting functions
scripts/Reanalysis_Exp1/    Joint Stengård, Experiment 1, and supporting analyses
scripts/Exp2/              Experiment 2 analyses
data/Exp1/                 Experiment 1 trials, demographics, and intermediate results
data/Exp2/                 Experiment 2 trials, reference simulations, and posterior draws
data/Stengard/             Stengård trials and intermediate results
data/sirota.csv            Supporting descriptive dataset
results/Exp2_PPC_r05/       Saved Experiment 2 predictive-check results and heatmaps
fig/png/, fig/tiff/         Supplied figures and figure output destinations
Table/                     LaTeX tables
docs/analysis_pipeline.md   Execution order, dependencies, and reproducibility details
```

## Analyses

The repository contains code for the following analyses:

- **Behavioural analyses:** Probability judgements, deviations from Bayesian predictions, heuristic matching, and effects of presentation format and time pressure.
- **Computational modelling:** Implementation and comparison of Bayesian, heuristic, sampling, and hybrid models, including the Mixture Heuristic-Anchored Bayesian Sampler (HABS).
- **Model comparison:** Approximate Bayesian computation with random forests (ABC-RF), followed by variational Bayesian model and family selection for Stengård and Experiment 1.
- **Parameter estimation:** Neural-network ABC for the shared `M_L` and `MH` analyses; rejection ABC for Experiment 2 `M_L` inference.
- **Model evaluation:** Classifier diagnostics and posterior predictive checks (PPC).
- **Figures and tables:** Scripts used to generate the results presented in the manuscript and supplementary materials.

## Reproducing the Analyses

The analyses were developed in R. Start a fresh session from the repository root, with dependencies installed as described in the [environment guide](docs/analysis_pipeline.md#environment-and-execution). For example, the first-block analysis consumes the supplied Experiment 2 trials and posterior draws:

```sh
Rscript --vanilla scripts/Exp2/Exp2_First_Block.R
```

For a detailed description of the analysis workflows, required inputs, execution order, and current reproducibility limitations, see [Analysis Pipeline](docs/analysis_pipeline.md).

Supplied intermediate results support downstream analyses. Raw reconstruction inputs for both experiments are absent, and end-to-end reproduction has not been verified. Existing numbered figures are preserved; future Experiment 2 exports use descriptive filenames. See the guide for [output mappings and numbering limitations](docs/analysis_pipeline.md#figures-and-tables).

## Citation

If you use this code or data, please cite the associated manuscript:

Lin, Y., Sanborn, A. N., & Chater, N. *Integrating Heuristics and Mental Sampling in Belief Updating*. Manuscript submitted to *Journal of Experimental Psychology: General*.

*Full publication details will be updated upon acceptance.*
