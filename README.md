# longtermLinkage

[![R-CMD-check](https://github.com/mathcg/longterm-linkage/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mathcg/longterm-linkage/actions/workflows/R-CMD-check.yaml)
[![DOI](https://img.shields.io/badge/DOI-10.1080%2F01621459.2025.2526703-blue)](https://doi.org/10.1080/01621459.2025.2526703)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

`longtermLinkage` implements survival-analysis methods for estimating long-term
treatment effects when randomized clinical-trial records are only partially
linked to observational follow-up data.

The package accompanies:

> Cheng, G., Chen, Y.-C., Unger, J. M., Till, C., and Zhao, Y.-Q. (2026).
> Long-Term Effect Estimation When Combining Clinical Trial and Observational
> Follow-Up Datasets. *Journal of the American Statistical Association*,
> 121(553), 492–501.
> [doi:10.1080/01621459.2025.2526703](https://doi.org/10.1080/01621459.2025.2526703)

## Methods

- Complete-case analysis (CC)
- Augmented complete-case analysis (CC+)
- Non-linked-as-censored analysis (NLAC)
- Inverse-probability-of-linkage weighting (IPLW)
- Augmented IPLW for the paper's supported simulation settings
- Cox models with time-independent or time-dependent treatment effects

## Installation

After the repository is published, install the development version with:

```r
# install.packages("remotes")
remotes::install_github("mathcg/longterm-linkage")
```

## Quick start

```r
library(longtermLinkage)

# Generate the first time-dependent simulation setting from the paper.
dat <- simulate_linkage_data(n = 1000, setting = 1, seed = 2026)

# Fit the proposed IPLW estimator.
fit <- fit_linkage_cox(
  dat,
  method = "iplw",
  cox_covariates = c("Z_1", "Z_2", "Z_3"),
  linkage_covariates = c("Z_1", "Z_2"),
  time_dependent = TRUE
)
fit
confint(fit)
```

## Synthetic PCPT example

The repository contains a synthetic dataset modeled after the linked Prostate
Cancer Prevention Trial analysis. It contains no original participant records.

```r
pcpt <- read_pcpt_synthetic()
pcpt_data <- prepare_pcpt_linkage_data(pcpt)
results <- compare_pcpt_methods(pcpt_data, time_dependent = TRUE)
results$beta_hat
```

See [`examples/`](examples/) for complete scripts. The exact files submitted
with the journal article are preserved under [`analysis/original/`](analysis/original/).

## Citation

If you use this software or these methods, please cite the JASA article:

```bibtex
@article{cheng2026longterm,
  author  = {Cheng, Gang and Chen, Yen-Chi and Unger, Joseph M. and Till, Cathee and Zhao, Ying-Qi},
  title   = {Long-Term Effect Estimation When Combining Clinical Trial and Observational Follow-Up Datasets},
  journal = {Journal of the American Statistical Association},
  volume  = {121},
  number  = {553},
  pages   = {492--501},
  year    = {2026},
  doi     = {10.1080/01621459.2025.2526703}
}
```

## Data note

The included PCPT dataset is synthetic and intended for software demonstration.
The original PCPT and Medicare-linked data are not distributed by this package.
