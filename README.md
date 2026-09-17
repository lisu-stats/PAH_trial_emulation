# Target Trial Emulation of Selexipag Initiation in Pulmonary Arterial Hypertension




## Overview

This repository contains the analysis code used to evaluate the long-term effect of Selexipag initiation on all-cause mortality among patients with idiopathic and heritable pulmonary arterial hypertension (PAH) 
using the target trial emulation framework.

The study compares two treatment strategies:

1. Initiate Selexipag when clinically indicated.
2. Never initiate Selexipag.

The primary analysis uses:

- Marginal structural models (MSMs) estimated by longitudinal targeted maximum likelihood estimation (LTMLE) with Super Learner ensemble for nuisance model estimation
- Delete-a-group jackknife for inference

---

## Associated Manuscript

**Long-term effects of Selexipag initiation in pulmonary arterial hypertension: a target trial emulation**

Manuscript currently in preparation.

---

## Repository Structure

```text
.
├── R/
│   ├── 0_PAH_datamanip.R
│   ├── 1_PAH_LTMLE_preparation.R
│   ├── 2_PAH_LTMLE_modelfit.R
│   ├── 3_PAH_group_jackknife.R
│   ├── 4_PAH_sensitivity_analysis.R
│   ├── 5_PAH_forest_plots.R
│   └── 6_PAH_risk_ratios.R
├── CITATION.cff
├── LICENSE
└── README.md
```

### Script descriptions

| Script | Purpose |
|----------|----------|
| `0_PAH_datamanip.R` | Data cleaning, harmonisation, and longitudinal dataset construction |
| `1_PAH_LTMLE_preparation.R` | Preparation of covariates, treatment indicators, censoring variables, and LTMLE inputs |
| `2_PAH_LTMLE_modelfit.R` | Main LTMLE analyses and marginal structural model estimation; subgroup analyses |
| `3_PAH_group_jackknife.R` | Delete-a-group jackknife variance estimation and confidence intervals; create results tables |
| `4_PAH_sensitivity_analysis.R` | Sensitivity analyses using E-values |
| `5_PAH_forest_plots.R` | Generation of manuscript and supplementary figures |
| `6_PAH_risk_ratios.R` | Calculation of risk ratios from primary analysis results |

---

## Methods

The primary analysis evaluates the effect of initiating Selexipag when clinically indicated compared with never initiating Selexipag.

Analyses were conducted using:

- Marginal structural models (MSMs) estimated by longitudinal targeted maximum likelihood estimation (LTMLE) with Super Learner ensemble for nuisance model estimation
- Delete-a-group jackknife for inference


Secondary analyses evaluated:

- Time-varying treatment effects
- Effect modification by GRIPHON trial eligibility
- Effect modification by incident vs prevalent case status
- Effect modification by Royal Papworth Hospital (RPH) status
- Sensitivity to unmeasured confounding using E-values

---

## Data Availability

The patient-level data used in this project cannot be shared publicly because of ethical, governance, and privacy restrictions.

This repository contains analysis code only.

Data are available upon reasonable request and subject to local laws, approvals and institutional requirements, for more information contact Professor Mark Toshner (mrt34@cam.ac.uk).

---

## Software Requirements

The analyses were conducted in R (4.6).

Main packages include:

```r
ltmle
SuperLearner
dplyr
ggplot2
parallel
```
Additional packages may be required depending on the local computing environment.


---

## Reproducibility

The analysis scripts are intended to be run sequentially:

```r
source("R/0_PAH_datamanip.R")
source("R/1_PAH_LTMLE_preparation.R")
source("R/2_PAH_LTMLE_modelfit.R")
source("R/3_PAH_group_jackknife.R")
source("R/4_PAH_sensitivity_analysis.R")
source("R/5_PAH_forest_plots.R")
source("R/6_PAH_risk_ratios.R")
```



---

## Citation

If you use this repository, please cite:

1. The associated manuscript.
2. This GitHub repository.

Citation metadata are provided in `CITATION.cff`.

---

## License

This project is licensed under the MIT License.

See the `LICENSE` file for details.

---

## Contact

Li Su

For questions regarding the analysis code, please open a GitHub Issue or contact the repository maintainer.
