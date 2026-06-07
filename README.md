# Adaptive-Ensemble-Learning-for-Short-Panel-Causal-Policy-Evaluation

This repository contains replication files for the paper:

**Adaptive Ensemble Learning for Short-Panel Causal Policy Evaluation**

The paper proposes perturbed cross-validation (PCV), an adaptive ensemble method for short-panel causal policy evaluation. PCV is designed for settings where the time dimension is short, the validation window is small, and researchers face multiple credible counterfactual estimators.

## Repository structure

- `data_raw/`: Original data files.
- `data_clean/`: Cleaned analysis-ready data.
- `code/`: R scripts for data cleaning, simulations, application, tables, and figures.
- `output/`: Generated tables, figures, simulation summaries, and application results.
- `paper/`: LaTeX source files and paper figures.
- `docs/`: Variable dictionary and replication notes.

## Replication

To reproduce the main tables and figures, run:

```r
source("master.R")
