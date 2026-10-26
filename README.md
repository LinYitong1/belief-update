# Integrating Heuristics and Mental Sampling in Belief Updating

This repository contains the data analysis and computational modelling code for the research project **Integrating Heuristics and Mental Sampling in Belief Updating**, submitted to *Journal of Experimental Psychology: General*.

## Overview

This project investigates how people combine simple heuristics and Bayesian sampling when updating beliefs in response to new evidence.

We examine whether probability judgements can be explained by a hybrid computational framework in which people use heuristics to generate initial estimates and refine them through limited Bayesian sampling.

The repository includes analyses from three datasets:

1. **Stengård dataset reanalysis:** Reanalysis of an existing dataset to evaluate heuristic, Bayesian, sampling, and hybrid models.
2. **Experiment 1:** An experimental investigation of belief updating using frequency and probability presentation formats.
3. **Experiment 2:** A follow-up experiment examining how time pressure influences probability judgements, heuristic use, and Bayesian sampling.

## Analyses

The repository contains code for the following analyses:

- **Behavioural analyses:** Probability judgements, deviations from Bayesian predictions, heuristic matching, and effects of presentation format and time pressure.
- **Computational modelling:** Implementation and comparison of Bayesian, heuristic, sampling, and hybrid models, including the Mixture Heuristic-Anchored Bayesian Sampler (HABS).
- **Model comparison:** Approximate Bayesian computation with random forests (ABC-RF), followed by variational Bayesian model and family selection for Stengård and Experiment 1.
- **Parameter estimation:** ABC for the shared `M_L` and `MH` analyses.
- **Model evaluation:** Classifier diagnostics and posterior predictive checks (PPC).
- **Figures and tables:** Scripts used to generate the results presented in the manuscript and supplementary materials.

## Reproducing the Analyses

The analyses were developed in R. Start a fresh session from the repository root, with dependencies installed as described in the [environment guide](docs/analysis_pipeline.md#environment-and-execution). For example, the first-block analysis consumes the supplied Experiment 2 trials and posterior draws:

```sh
Rscript --vanilla scripts/Exp2/Exp2_First_Block.R
```

## Citation

If you use this code or data, please cite the associated manuscript:

Lin, Y., Sanborn, A. N., & Chater, N. *Integrating Heuristics and Mental Sampling in Belief Updating*. Manuscript submitted to *Journal of Experimental Psychology: General*.

*Full publication details will be updated upon acceptance.*
