\name{fit_linkage_cox}
\alias{fit_linkage_cox}
\title{Fit a Cox Model Under a Linkage-Analysis Strategy}
\usage{
fit_linkage_cox(
  data,
  method = c("iplw", "nlac", "complete_case", "augmented_complete_case", "oracle"),
  cox_covariates,
  linkage_covariates = NULL,
  time_dependent = FALSE
)
}
\arguments{
\item{data}{A list containing a subject-level \code{df} data frame and, for a
time-dependent analysis, a counting-process \code{tdata} data frame.}
\item{method}{One of \code{"oracle"}, \code{"complete_case"},
\code{"augmented_complete_case"}, \code{"nlac"}, or \code{"iplw"}.}
\item{cox_covariates}{Character vector naming Cox-model covariates.}
\item{linkage_covariates}{Character vector naming linkage-model covariates;
required for \code{method = "iplw"}.}
\item{time_dependent}{Whether to fit a counting-process Cox model.}
}
\value{
An object of class \code{longterm_fit} with coefficients, standard errors,
covariance matrix, and the underlying model result.
}
\description{
Fit the oracle, complete-case (CC), augmented complete-case (CC+),
non-linked-as-censored (NLAC), or inverse-probability-of-linkage weighted
(IPLW) analysis described by Cheng et al. (2026).
}
\examples{
\donttest{
dat <- simulate_linkage_data(n = 500, setting = 1, seed = 2026)
fit <- fit_linkage_cox(
  dat,
  method = "iplw",
  cox_covariates = c("Z_1", "Z_2", "Z_3"),
  linkage_covariates = c("Z_1", "Z_2"),
  time_dependent = TRUE
)
fit
}
}
