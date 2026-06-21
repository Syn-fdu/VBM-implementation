# Variance-Based Sensitivity Analysis for Hidden Confounding

[中文说明](README.zh-CN.md) | English


The guiding question is:

> How strong would unobserved confounding have to be to overturn the conclusion that high fish/shellfish consumption is associated with higher blood mercury?


---

## What is included


1. loads the NHANES fish-consumption data;
2. estimates the ATT using inverse probability weighting;
3. computes VBM sensitivity thresholds and covariate benchmarks;
4. compares VBM with marginal sensitivity model (MSM) benchmarks;
5. automatically generates standard and enhanced figures.

The four extension scripts are:

| Script | Purpose |
|---|---|
| `extension_01_hidden_strength_vbm_msm.R` | Compares VBM and MSM under different hidden-confounder strengths. |
| `extension_02_vbm_ps_misspecification.R` | Studies VBM under observed propensity-score model misspecification. |
| `extension_03_good_overlap_instability.R` | Stress-tests denominator instability under very good overlap. |
| `extension_04_rgm_vbm_msm.R` | Explores an RGM / TV-L1 sensitivity model and compares it with VBM and MSM. |

---

## Repository structure

```text
.
├── setup.R
├── extension_01_hidden_strength_vbm_msm.R
├── extension_02_vbm_ps_misspecification.R
├── extension_03_good_overlap_instability.R
├── extension_04_rgm_vbm_msm.R
├── functions/
├── rgm_model/
├── scripts/
└── output/
```

There is intentionally **no `extensions/` source folder**. The four extension scripts are kept at the repository root so that readers can find and run them quickly.

---

## Quick start

Clone the repository, open the project folder in R, and run:

```r
source("setup.R")   # optional: install missing R packages
```


3. all four extension scripts;
4. each extension's own enhanced plotting script.

Outputs are regenerated under:

```text
output/tables/
output/figures/
output/figures_enhanced/
output/extension_results/
```


---

## Run one extension only

```r
source("extension_01_hidden_strength_vbm_msm.R")
source("extension_02_vbm_ps_misspecification.R")
source("extension_03_good_overlap_instability.R")
source("extension_04_rgm_vbm_msm.R")
```

---

## Configuration

Settings are in:

```text
functions/config.R
```

Important defaults:

```r
run_extensions = TRUE
run_extension_01_hidden_strength_vbm_msm = TRUE
run_extension_02_vbm_ps_misspecification = TRUE
run_extension_03_good_overlap_instability = TRUE
run_extension_04_rgm_vbm_msm = TRUE
```

Set any switch to `FALSE` for a shorter run.

---

## Requirements

R packages:

- `CrossScreening`
- `survey`
- `boot`
- `ggplot2`
- `dplyr`
- `tidyr`
- `readr`
- `scales`

Python packages for enhanced plots:

- `pandas`
- `numpy`
- `matplotlib`

Install Python packages with:

```bash
pip install pandas numpy matplotlib
```

---


Sensitivity analysis makes hidden assumptions visible. This project compares three ways to describe hidden confounding:

- **VBM / L2-R2**: measures systematic residual weight variation.
- **RGM / L1-TV**: measures how much probability mass hidden confounding can move.


---

## Reference

Huang & Pimentel (2025). Variance-based sensitivity analysis. *Biometrika*.
