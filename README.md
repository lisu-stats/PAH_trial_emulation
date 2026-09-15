# Target Trial Emulation of Selexipag Initiation in Pulmonary Arterial Hypertension

![License](https://img.shields.io/badge/license-MIT-greeng.shields.io/badge/language-R-blue

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
│   ├── 01_data_preparation.R
│   ├── 02_primary_analysis.R
│   ├── 03_secondary_analyses.R
│   ├── 04_evalue_analysis.R
│   └── helper_functions.R
├── output/
│   ├── figures/
│   └── tables/
├── README.md
├── LICENSE
└── CITATION.cff
```

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

Researchers interested in accessing the underlying data should contact Professor Mark Toshner (mrt34@cam.ac.uk) and comply with all applicable governance and ethics requirements.

---

## Software Requirements

The analyses were conducted in R.

Main packages include:

```r
ltmle
SuperLearner
data.table
dplyr
tidyverse
ggplot2
parallel
```

---

## Reproducibility

The scripts are intended to be run sequentially:

```r
source("R/01_data_preparation.R")
source("R/02_primary_analysis.R")
source("R/03_secondary_analyses.R")
source("R/04_evalue_analysis.R")
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
