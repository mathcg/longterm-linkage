\name{run_linkage_simulation}
\alias{run_linkage_simulation}
\title{Run One Simulation Replicate from the Paper}
\usage{
run_linkage_simulation(
  n = 2000L,
  seed = 1L,
  method = c(
    "iplw", "nlac", "complete_case", "augmented_complete_case", "oracle", "aiplw"
  ),
  setting = 1L,
  time_dependent = setting <= 5L,
  n_boot = 100L
)
}
\arguments{
\item{n}{Positive integer sample size.}
\item{seed}{Positive integer random seed.}
\item{method}{One of \code{"oracle"}, \code{"complete_case"},
\code{"augmented_complete_case"}, \code{"nlac"}, \code{"iplw"}, or
\code{"aiplw"}.}
\item{setting}{Integer from 1 through 7 selecting a paper simulation design.}
\item{time_dependent}{Whether to use the time-dependent analysis.}
\item{n_boot}{Number of bootstrap samples for AIPLW.}
}
\value{
A list with coefficient estimates, standard errors, and coverage indicators
for the known simulation truth.
}
\description{
Run one estimator on one generated dataset using the simulation designs from
the main paper and supplementary material.
}
\examples{
\donttest{
result <- run_linkage_simulation(
  n = 500, seed = 2026, method = "nlac", setting = 1
)
result$beta
}
}
