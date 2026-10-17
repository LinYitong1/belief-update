# Analysis pipeline and file map

This guide describes the current source and supplied files on 17 October 2026. Execution order follows code readers, writers, and helper calls; it is not a verified history of the saved numerical results. The consistency cleanup changes paths, an explicit helper import, figure export names, and documentation only. No analyses or outputs were regenerated and no data or existing result files were moved, renamed, or deleted during this cleanup.

All paths are relative to the repository root. Bare data filenames in the shared workflow refer to `data/Exp1/` for Experiment 1 and `data/Stengard/` for Stengård; bare data filenames in the Experiment 2 workflow refer to `data/Exp2/`. Supporting `data/sirota.csv` remains at its existing location.

## Study mapping and glossary

| Study / supporting dataset | Input and scripts | Scope |
| --- | --- | --- |
| Stengård reanalysis | `data/Stengard/s.csv`; `scripts/Reanalysis_Exp1/` | Presentation-format comparisons and computational model comparison/inference. Code uses `stengard`; external-data provenance and full citation need confirmation. |
| Experiment 1 | `data/Exp1/df.csv`; `scripts/Reanalysis_Exp1/` | Frequency/probability presentation formats and repeated judgments. Supplied trials are already cleaned. |
| Experiment 2: time pressure | `data/Exp2/Exp_2.csv`; `scripts/Exp2/` | Self-paced/time-limited blocks, first-block and two-block analyses, and block order. |
| Sirota supporting descriptive dataset | `data/sirota.csv`; shared `00_plot.R` and `00_stats_table.R` | Figure 1 identifies Sirota et al. (2024). The table script also loads this input, although its exported tables use the other datasets. Full citation and preparation history need confirmation. |

See [shared plots](../scripts/Reanalysis_Exp1/00_plot.R), [dataset loaders](../R/plot_helpers.R), and [Experiment 2 cleaning](../scripts/Exp2/Exp2_Data_Clean.R).

| Term | Implemented meaning |
| --- | --- |
| `M_L` | `large_mixed_model`: mixtures of heuristic-anchored Bayesian sampling, base heuristic responses, and random responses; six heuristic weights. Display names vary across plotting and model-comparison helpers. |
| `MH` | `mixed_heuristic_model`: trial-wise selection among six heuristics and a random-response component. |
| ABC-RF | Approximate Bayesian computation with random forests for comparison of 20 candidate models in the shared workflow. This stage is followed by variational Bayesian model/family selection. |
| ABC parameter estimation | The shared posterior helper uses `abc::abc(method = "neuralnet")`; Experiment 2 uses rejection ABC. Both use tolerance 0.05. |
| PPC | Posterior predictive check: simulate responses from estimated posterior parameters and compare response statistics with observations. |
| PPP | The implemented two-sided predictive p-value is `1 - 2 * abs(u - 0.5)`, where `u` is the proportion of simulated statistics at least as large as the observation. It is not a model posterior probability. |

Definitions follow [models](../R/models.R), [candidate model calls](../scripts/Reanalysis_Exp1/02_generate_newprior.R), [model comparison](../scripts/Reanalysis_Exp1/04_abcrf_individual_inference.R), [posterior helpers](../R/posterior_helper.R), and [Experiment 2 inference](../scripts/Exp2/Exp2_Post_Parameters.R). `_s`/`dts` designate Stengård objects. Preserve literal filenames such as `_r05` and `newprior`; a complete historical naming convention is not established.

## Reproduction status

No workflow has been reproduced end to end in this cleanup. “Inputs present” refers to file availability, not runtime success, package availability, or compatibility of cached objects.

| Workflow | Available starting point | Missing inputs / qualification |
| --- | --- | --- |
| Experiment 1 raw reconstruction | `data/Exp1/demographic.csv`; cleaned `df.csv` is supplied | `data/Exp1/data_raw/` is absent. `data/Exp1/renamed_data/` is generated staging, also absent. |
| Shared modelling rebuild | Both cleaned study CSVs; prior, prediction, human/simulation summary, and model-probability caches | All statically identified downstream file inputs are present; runtime unverified. |
| Shared posterior/PPC branches | Both trial CSVs, prior tables, simulation summaries, Stengård human summary, and model-probability caches | Scripts recompute inference and PPC; subject-result RDS caches are not supplied and are not explicitly exported by current code. |
| Shared plots and tables | Three analysis CSVs, two prediction caches, and four PPP PNGs | Inputs present. Plotting includes stochastic simulations; it is not merely an image-assembly step. |
| Experiment 2 raw cleaning | Supplied cleaned `data/Exp2/Exp_2.csv` supports a separate downstream route | Expected raw directory `data/Exp2/Exp_2/` is absent. Optional `data/Exp2/Pilot_demographic.csv` is absent. |
| Experiment 2 human/reference summaries | Cleaned trials | Inputs and cached outputs present. |
| Experiment 2 posterior estimation and PPC | Reference table, human summary, posterior draws, and cleaned trials as appropriate | Required files present. Human-summary spelling is now consistently `Exp2_Human_Stats.rds`. |
| Experiment 2 first-/two-block analyses | Cleaned trials and posterior draws | Inputs present. Explicit plotting helper import corrected; runtime unverified. |

No missing input has been replaced or fabricated. Raw-input paths now sit within the intended study directories, but no raw directories were created. External dataset preprocessing is not reconstructable from the maintained source alone.

## Environment and execution

Use the repository root as the working directory. Many scripts use literal `source("R/...")` and relative paths; others use `here::here()`. Open `belief-update_copy.Rproj` or start `Rscript --vanilla` from the root. The project retains default workspace/history settings; a fresh command-line session avoids implicit `.RData` restoration.

R 4.3.3 parses the source; this does not establish a validated runtime environment. Syntax uses features introduced in R 4.1. No package lock or complete tested installation recipe is supplied.

| Purpose | Declared packages |
| --- | --- |
| Data handling | `data.table`, `tidyverse`, `readxl`, `jsonlite` |
| Simulation and inference | `MCMCpack`, `abc`, `abcrf`, `ranger`, `matrixStats` |
| Statistical analyses | `lme4`, `lmerTest`, `emmeans`, `clubSandwich`, `afex`, `brms`, `binom` |
| Figures and tables | `ggplot2`, `patchwork`, `cowplot`, `ggpubr`, `gridExtra`, `ggforce`, `magick`, `knitr`, `kableExtra`, `gt`, `broom`, `broom.mixed`, `PairedData` |
| Execution support | `pacman`, `here`, `foreach`, `doParallel`, `future`, `furrr`, `progressr` |

Additional declarations include `readr`, `dplyr`, `tidyr`, `purrr`, `stringr`, `tibble`, `magrittr`, `forcats`, `scales`, `viridisLite`, `RColorBrewer`, `devtools`, and base/recommended packages such as `parallel`, `stats`, and `grid`; namespace references include `rlang`. These are source declarations, not evidence that every package feature is used or installed.

`pacman::p_load()` and the posterior scripts' package loader can install missing packages. The posterior loader calls `purrr::walk()` before ensuring its dependencies. Preinstall packages in the intended environment. Keep helper source order: loading helpers together can overwrite definitions. Separate fresh sessions are preferable. Scripts use fixed output paths and can overwrite matching results; reproduce in a separate copy when ready.

For example, this command consumes supplied Experiment 2 trials and posterior draws:

```sh
Rscript --vanilla scripts/Exp2/Exp2_First_Block.R
```

This invocation has not been executed as part of the cleanup. Script numbering is not a run-all order.

## Experiment 1 and Stengård reanalysis

All scripts below are in `scripts/Reanalysis_Exp1/`. They jointly process the two studies, so their reference and parameter inputs must remain matched.

```text
Exp1 raw inputs (absent) → 00_data_process.R → data/Exp1/df.csv
                                              + data/Stengard/s.csv
  → 01_simulate_parameter_newprior.R → 02_generate_newprior.R
  → 03_summary_stats_newprior.R → 04_abcrf_individual_inference.R
  → 05_posterior.R / 05_posterior_MH.R (independent branches)
  → PPP images consumed by 00_plot.R for Appendix F
```

| Script | Inputs | Outputs / downstream role |
| --- | --- | --- |
| [00_data_process.R](../scripts/Reanalysis_Exp1/00_data_process.R) | Missing `data/Exp1/data_raw/`; demographic summaries use `data/Exp1/demographic.csv` | Copies raw files into `data/Exp1/renamed_data/`; writes `data/Exp1/df.csv`; console summaries. |
| [01_simulate_parameter_newprior.R](../scripts/Reanalysis_Exp1/01_simulate_parameter_newprior.R) | Both study CSVs; selects primary subject `1` and Stengård participant-format `1_1` as designs | `data/Exp1/parameter_dt_r05.rds`, `data/Stengard/parameter_dt_s_r05.rds`. |
| [02_generate_newprior.R](../scripts/Reanalysis_Exp1/02_generate_newprior.R) | Both parameter tables; `R/models.R` | `prediction_dt_r05.rds` and `prediction_dt_s_r05.rds` in the respective study directories; 20 models, predictions rounded to two decimals. |
| [03_summary_stats_newprior.R](../scripts/Reanalysis_Exp1/03_summary_stats_newprior.R) | Both predictions and trial CSVs; `R/metrics.R` | Experiment 1 `Human_Summary_dt.rds` / `Simulate_Summary_dt_r05.rds`; Stengård `Human_Summary_dts.rds` / `Simulate_Summary_dts_r05.rds`; `fig/variance.png`. |
| [04_abcrf_individual_inference.R](../scripts/Reanalysis_Exp1/04_abcrf_individual_inference.R) | Both human/simulation summary pairs and trial CSVs; VB-BMS, ABC-RF, plotting helpers | `data/Exp1/posterior_long.rds`, `data/Stengard/posterior_long_s.rds`; Figures 13 and 6; model/family tables printed rather than exported. |
| [05_posterior.R](../scripts/Reanalysis_Exp1/05_posterior.R) | Both trial CSVs, both prior tables, both simulation summaries, Stengård human summary, both model-probability files | Recomputes `M_L` posteriors/PPC. Figure 8, descriptive Experiment 1 predictive-distribution exports, Appendices E2/E4, and PPP histograms. Result objects remain in memory. |
| [05_posterior_MH.R](../scripts/Reanalysis_Exp1/05_posterior_MH.R) | Same file set as the `M_L` branch; model assignments also support a Stengård subgroup analysis | Recomputes `MH` posteriors/PPC. Figure 7, Appendices E1/E3, PPP histograms, and `Table/slope_s_mh_non.tex`. Experiment 1 overlay remains in memory. |
| [00_plot.R](../scripts/Reanalysis_Exp1/00_plot.R) | Both study CSVs, `data/sirota.csv`, both prediction tables, four saved PPP PNGs | Figures 1–5 and 9–12, Appendices D/F; descriptive tests and stochastic model illustrations. |
| [00_stats_table.R](../scripts/Reanalysis_Exp1/00_stats_table.R) | Both study CSVs and `data/sirota.csv`; `R/table_helper.R` | `Table/Table2.tex`, `Table/Table4.tex`, `Table/Table5.tex`; additional overlap tables/summaries printed or retained in memory. |

The table script is independent of modelling caches. The plot script's `00_` prefix does not mean it should run first: it loads predictions at startup and later requires both posterior branches' PPP images. Supplied prediction/PPP files permit a cached entry route. Sirota Figure 1 has no separate standalone script.

Both posterior branches need all listed loaded inputs, including objects not subsequently used in estimation. They recompute inference/PPC; neither reads a saved `*subject_results*.rds` object. Do not regenerate priors independently of their matching reference simulations: row alignment is scientifically significant.

## Experiment 2

All scripts below are in `scripts/Exp2/`. Data and intermediate RDS destinations are under `data/Exp2/`.

```text
Exp_2/ raw sessions (absent) → Exp2_Data_Clean.R → Exp_2.csv
Exp_2.csv → Exp2_Human_Stats.R ───────┐
Exp_2.csv → Exp2_Prior_Parameters.R ──┴→ Exp2_Post_Parameters.R
                                      → posterior draws
posterior draws + cleaned trials → first-block / two-block analyses
posterior draws + trials + reference + human summaries → PPC
```

| Script | Inputs | Outputs / downstream role |
| --- | --- | --- |
| [Exp2_Data_Clean.R](../scripts/Exp2/Exp2_Data_Clean.R) | Missing `data/Exp2/Exp_2/` session directory; optional `Pilot_demographic.csv` | `Exp_2.csv`; completion/catch-trial and trial-validity checks, timing reconstruction, console diagnostics. |
| [Exp2_Human_Stats.R](../scripts/Exp2/Exp2_Human_Stats.R) | Cleaned trials; metrics and mean–variance helpers | `Exp2_Human_Stats.rds`; both blocks required, no minimum 45-trial filter. |
| [Exp2_Prior_Parameters.R](../scripts/Exp2/Exp2_Prior_Parameters.R) | Cleaned trials; models, metrics, mean–variance helpers | First complete 54-trial block (27 items repeated twice) supplies the template; 10,000 `M_L` draws; `Exp2_ML_Parameters_r05.rds` and `Exp2_ML_reference_r05.rds`. |
| [Exp2_Post_Parameters.R](../scripts/Exp2/Exp2_Post_Parameters.R) | Reference table and human summaries | Rejection ABC for 11 parameters; `Exp2_ML_Posterior_Draws_r05.rds`; paired tests and summaries printed. |
| [Exp2_PPC.R](../scripts/Exp2/Exp2_PPC.R) | Cleaned trials, posterior draws, reference table, human summaries; models/metrics/mean–variance helpers | `results/Exp2_PPC_r05/Exp2_PPC_Results.rds`; full- and first-block console summaries. Existing heatmap PNGs have no explicit writer in the current script. |
| [Exp2_First_Block.R](../scripts/Exp2/Exp2_First_Block.R) | Cleaned trials and posterior draws | First-block behavioural and parameter analyses; console/in-memory output, no explicit export. |
| [Exp2_Two_Blocks.R](../scripts/Exp2/Exp2_Two_Blocks.R) | Cleaned trials and posterior draws; metrics, models, posterior and plotting helpers | Two-block/order analyses; descriptive Experiment 2 behaviour and parameter figure pairs listed below; console/in-memory results. |

Human-summary and simulated-reference generation are independent after cleaned trials exist. The three downstream branches are independent of one another. PPC additionally needs the human summary and reference. All references now use the supplied uppercase `Exp2_Human_Stats.rds`; the plotting import now resolves to `R/plot_helpers.R`.

## Figures and tables

Figure exports use `fig/png/` and `fig/tiff/`, except the existing summary-stage diagnostic destination `fig/variance.png`. No supplied image or table was changed.

| Producer | Code-defined exports | Availability / qualification |
| --- | --- | --- |
| Shared `00_plot.R` | Figures 1–5, 9–12; Appendices D/F | PNGs present. TIFFs present except `Figure_5.tiff`; `Figure_5_r05.tiff` exists, equivalence unverified. |
| Shared `04_abcrf_individual_inference.R` | Figures 6 and 13 | Both formats present. |
| Shared `05_posterior.R` | Figure 8; `Exp1_ML_Predictive_Distributions`; Appendices E2/E4; `PPP_hist_Stengard_indiv` and `PPP_hist_Experimental_indiv` | Numbered/appendix/PPP files present. Descriptive Experiment 1 pair is a new export destination and has not been generated. |
| Shared `05_posterior_MH.R` | Figure 7; Appendices E1/E3; `PPP_hist_Stengard_indiv_MH` and `PPP_hist_exp_indiv_MH` | Both formats present. Existing Figure 14 has no explicit writer in the current script, which builds the Experiment 1 overlay in memory. |
| `Exp2_Two_Blocks.R` | `fig/png/Exp2_Behaviour_Block_Order.png`, `fig/tiff/Exp2_Behaviour_Block_Order.tiff`, `fig/png/Exp2_Parameters_Block_Order_r05.png`, `fig/tiff/Exp2_Parameters_Block_Order_r05.tiff` | Descriptive exports have not been generated. Supplied Figures 15/16 remain untouched. |
| Shared `00_stats_table.R` | `Table/Table2.tex`, `Table/Table4.tex`, `Table/Table5.tex` | All three present. Tables 1 and 3 are described in source as manually created in LaTeX. |
| Shared `05_posterior_MH.R` | `Table/slope_s_mh_non.tex` | Present. Other posterior/overlap summaries are not all exported. |
| `Exp2_PPC.R` | `results/Exp2_PPC_r05/Exp2_PPC_Results.rds` | Present; `PPC_Full_Heatmap.png` and `PPC_First_Block_Heatmap.png` also exist, but current source does not export them. |

Visual inspection identifies supplied `fig/png/Figure_15.png` as Experiment 2 behaviour and `Figure_16.png` as Experiment 2 parameters. Previously, the Experiment 1 `M_L` script targeted `Figure_15`, which would overwrite unrelated content. Its destination is now `Exp1_ML_Predictive_Distributions`; Experiment 2 uses the descriptive names above. No manuscript figure was renumbered. Exact manuscript correspondence and historical TIFF/PNG provenance still require author confirmation; existing numbered assets alone do not validate their generating code or numerical equivalence.

`fig/variance.png` is absent. `Apendix_D.tiff` / `Appendix_D.tiff` and the corresponding Appendix F spellings are byte-identical duplicate pairs; all were preserved.

## Shared helper inventory

Helpers are sourced by entry scripts. Preserve source order and calling context.

| File | Responsibility |
| --- | --- |
| [models.R](../R/models.R) | Bayesian, heuristic, sampling, and mixture response models. |
| [metrics.R](../R/metrics.R) | Matching/deviation statistics, regression summaries, and repeated-response variance. |
| [mean_variance.R](../R/mean_variance.R) | Repeated-item mean–variance slope. |
| [abcrf_helpers.R](../R/abcrf_helpers.R) | Custom random-forest leaf/probability calculations. |
| [VB_bms.R](../R/VB_bms.R) | Variational Bayesian model/family selection and reporting. |
| [posterior_helper.R](../R/posterior_helper.R) | ABC inference, predictive simulation, PPP and parameter summaries/plots. |
| [plot_helpers.R](../R/plot_helpers.R) | Dataset loaders and plotting, including Experiment 2 exports. Some helpers depend on caller/global objects. |
| [table_helper.R](../R/table_helper.R) | Dataset loaders, statistical summaries, and LaTeX exports. Its Stengård identifier handling differs from the plotting loader. |

Dataset loaders accept paths from callers; their transformation logic was not changed. The Experiment 2 figure-saving helper now targets existing format subdirectories. Other shared helpers required no edits.

## Storage roles and historical artifacts

| Location | Current contents |
| --- | --- |
| `data/Exp1/` | `df.csv`, `demographic.csv`, and five RDS caches: parameters, predictions, human/simulation summaries, model probabilities. |
| `data/Stengard/` | `s.csv` and five corresponding RDS caches with Stengård suffixes. |
| `data/Exp2/` | `Exp_2.csv`, `Exp2_Human_Stats.rds`, `Exp2_ML_Parameters_r05.rds`, `Exp2_ML_reference_r05.rds`, `Exp2_ML_Posterior_Draws_r05.rds`. |
| `data/sirota.csv` | Supporting descriptive input; remains at the root of `data/`. |
| `results/Exp2_PPC_r05/` | PPC result RDS and two supplied heatmaps. |
| `fig/`, `Table/` | Supplied presentation artifacts; current producer/availability qualifications above. |
| `.RData`, `.Rhistory`, `.Rproj.user/` | Local session/recovery material; not an environment specification or required analysis input. |

The previously documented eight `boardprior` caches and four `*subject_results*.rds` caches are absent from this snapshot. No replacements were created. Neither shared posterior script explicitly saves subject-result objects. The formerly descriptive Experiment 2 figure files are also absent; current exports now use descriptive names without altering supplied numbered images.

## Randomness and computational requirements

| Stage | Existing controls / scale |
| --- | --- |
| Shared priors | Seed 123 within `run_simulations`; 10,000 draws per design; forked `mclapply`, detected cores minus two, minimum one. |
| Shared predictions | Stochastic calls across 20 models; no top-level seed. |
| RF/model selection | 500 trees; no top-level seed; parallel backend uses detected cores minus four, minimum one; BMS exceedance calculations request 1,000,000 samples. |
| Shared `M_L` posterior/PPC | Top-level seed 20261011; neural-network ABC, tolerance 0.05. |
| Shared `MH` posterior/PPC | No top-level seed; helper-specific seeds include plotting jitter. |
| Experiment 2 references | Seed 123; 10,000 draws; forked multicore execution. |
| Experiment 2 posterior | Rejection ABC, tolerance 0.05; iteration alignment is checked. |
| Experiment 2 PPC | Seed 20261011; 500 posterior draws per participant/block; 41 reference statistics expected by this design. |
| Descriptive plots | Include stochastic simulations without a top-level seed. |

Seeds do not fix parallel RNG streams, worker counts, package algorithms, or historical interactive-session state. Cross-platform/bitwise reproduction has not been tested. Forked execution needs consideration on Windows. Runtime and memory have not been benchmarked; reference tables, posterior objects, and TIFF exports can be large. Do not run a blanket run-all loop without checking resource use and fixed output destinations.

## Unresolved interpretation and reproducibility issues

These qualifications describe existing implementation; the cleanup does not change scientific assumptions.

- Both experiments' raw-input directories are absent; supplied cleaned trials permit downstream entry. External dataset preparation/provenance, full citations, redistribution terms, and public-data de-identification require confirmation.
- The two `load_clean_data()` definitions handle Stengård identifiers differently. Sourcing the table helper after the plotting helper changes the active definition. They were not consolidated or reordered.
- RF quantities labelled OOB include predictions on supplied training data; they should not automatically be interpreted as held-out recovery validation.
- Experiment 2 human summaries require two blocks but do not impose a 45-trial minimum. No filter was added.
- `Exp2_Post_Parameters.R` prints a BH-correction label for paired tests but does not apply `p.adjust` in that block. No correction or inferential method was changed.
- Experiment 2 human-summary and PPC code use different minimum repeated-item counts for mean–variance slopes. Shared reference predictions are rounded while shared posterior-predictive simulations are not uniformly rounded. These conventions were retained.
- Some results remain in memory or console output. There is no complete scripted report export, maintained R Markdown rendering workflow, or automated regression suite identified.
- Figure numbering/provenance qualifications are recorded above. No numerical equivalence between supplied figures and current source is asserted.
- No repository licence or locked runtime environment is supplied. Local recovery artifacts and demographic/identifier fields remain in the working copy; this cleanup does not establish release permissions or anonymisation.

Before any analytical refactoring, preserve matched datasets, priors/reference tables, assignments, posterior draws, and checks; agree on deterministic and stochastic comparison criteria. Existing results have not been regenerated or certified.

## Validation history

Earlier documentation recorded parse/source and hash checks for an older layout. Those historical file counts and availability claims do not describe this snapshot.

For the current consistency cleanup, validation is static: parse every maintained R file without sourcing it; resolve explicit helper imports; inspect literal and constructed file paths and output destinations; verify documentation links; compare source edits and hashes against a pre-cleanup snapshot. Raw-input and optional-demographic exceptions are listed in the reproduction-status table. The model implementations, statistical expressions, priors, inference methods, and numerical assets are outside the permitted changes.

No package installation, full pipeline, simulation, model fit, plot, table, or independent review was executed as part of this cleanup. Runtime and numerical reproduction remain unverified.

Current check results: all 24 R files parse with R 4.3.3; all 31 explicit `source()` paths resolve with exact filename case. Static tracing checked 62 resolved input references, including literal paths, `here::here()`, `file.path()`, and caller-supplied loader paths. Only the raw-input/staging directories and optional demographic file listed above are absent. Parameterised helper readers were traced to their callers; session CSV contents cannot be checked without the absent raw directories. All 36 relative documentation links/anchors resolve. Data outputs route to their study folders; PPC results, figures, and tables retain their dedicated output roots.

Comparison of parsed expressions before/after, normalising only the authorised path, import, export-name, and saved-file-message substitutions, found no other code changes. SHA-256 comparison confirmed all supplied data, results, figures, and tables are byte-identical to the pre-cleanup snapshot. The cleanup modified 16 analysis scripts, one plotting helper, and two documentation files; no non-session files were added or deleted. Live `.Rproj.user/` state changed separately during the checks and was excluded from the final non-session manifest comparison; it was not edited by the cleanup.
