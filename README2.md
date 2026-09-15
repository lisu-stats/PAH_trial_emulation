# Target Trial Emulation of Selexipag Initiation in Pulmonary Arterial Hypertension

## Overview

This repository contains the analysis code used to evaluate the long-term effect of Selexipag initiation on all-cause mortality among patients with idiopathic and heritable pulmonary arterial hypertension (PAH) using a target trial emulation framework.

The study compares two treatment strategies:

1. Initiate Selexipag when clinically indicated.
2. Never initiate Selexipag.

The primary analysis uses:

- Longitudinal Targeted Maximum Likelihood Estimation (LTMLE)
- Super Learner ensemble modelling
- Marginal Structural Models (MSMs)
- Delete-a-group jackknife inference

---

## Associated Manuscript

**Long-term effects of Selexipag initiation in pulmonary arterial hypertension: a target trial emulation**

*Manuscript currently under review.*

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
│   ├── 5_make_forest_plots.R
│   └── 6_PAH_risk_ratios.R
├── CITATION.cff
├── LICENSE
└── README.md
```

---

## Methods

The primary target trial emulation evaluates the effect of initiating Selexipag when clinically indicated compared with never initiating Selexipag.

Analyses were conducted using:

- Longitudinal Targeted Maximum Likelihood Estimation (LTMLE)
- Super Learner
- Marginal Structural Models (MSMs)
- Grouped (delete-a-group) jackknife variance estimation

Secondary analyses evaluated:

- Treatment-effect heterogeneity by follow-up duration
- Effect modification by GRIPHON trial eligibility
- Effect modification by incident versus prevalent case status
- Effect modification by Royal Papworth Hospital (RPH) status
- Sensitivity to unmeasured confounding using E-values

---

## Workflow

The analysis scripts are intended to be run sequentially:

```r
source("R/0_PAH_datamanip.R")
source("R/1_PAH_LTMLE_preparation.R")
source("R/2_PAH_LTMLE_modelfit.R")
source("R/3_PAH_group_jackknife.R")
source("R/4_PAH_sensitivity_analysis.R")
source("R/5_make_forest_plots.R")
source("R/6_PAH_risk_ratios.R")
```

### Script descriptions

#### 0_PAH_datamanip.R
Data cleaning, harmonisation, and construction of the longitudinal analysis dataset.

#### 1_PAH_LTMLE_preparation.R
Preparation of treatment, censoring, outcome, and time-varying covariate histories required for LTMLE.

#### 2_PAH_LTMLE_modelfit.R
Primary LTMLE analysis and marginal structural model estimation.

#### 3_PAH_group_jackknife.R
Delete-a-group jackknife variance estimation and confidence interval construction.

#### 4_PAH_sensitivity_analysis.R
Secondary analyses including:
- Time-varying treatment effects
- GRIPHON eligibility subgroup analyses
- Incident versus prevalent case analyses
- Royal Papworth Hospital sensitivity analyses

#### 5_make_forest_plots.R
Generation of manuscript and supplementary figures.

#### 6_PAH_risk_ratios.R
Calculation of risk ratios and E-values from estimated risks and risk differences.

---

## Data Availability

Patient-level data are not publicly available because of ethical, governance, and privacy restrictions.

This repository contains analysis code only.

Researchers interested in accessing the underlying data should contact the study investigators and comply with all relevant governance and ethics requirements.

---

## Software Requirements

The analyses were conducted in R.

Major packages include:

```r
ltmle
SuperLearner
data.table
dplyr
ggplot2
parallel
```

Additional packages may be required depending on the local computing environment.

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

**Li Su**

For questions regarding the analysis code, please open a GitHub Issue or contact the repository maintainer.
