\name{simulate_linkage_data}
\alias{simulate_linkage_data}
\title{Simulate Linked Trial and Observational Follow-Up Data}
\usage{
simulate_linkage_data(n, setting = 1L, seed = 1L)
}
\arguments{
\item{n}{Positive integer sample size.}
\item{setting}{Integer from 1 through 7 selecting a paper simulation design.
Settings 1--5 use a time-dependent treatment effect; settings 6--7 use
time-independent covariates.}
\item{seed}{Positive integer random seed.}
}
\value{
A list with \code{df}, the subject-level data, and (for settings 1--5)
\code{tdata}, a counting-process representation.
}
\description{
Generate one dataset from a simulation design used in Cheng et al. (2026).
}
\examples{
dat <- simulate_linkage_data(n = 100, setting = 1, seed = 2026)
names(dat)
}
