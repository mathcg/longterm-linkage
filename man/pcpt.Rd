\name{pcpt_workflow}
\alias{read_pcpt_synthetic}
\alias{prepare_pcpt_linkage_data}
\alias{compare_pcpt_methods}
\title{Synthetic PCPT Linkage Analysis Workflow}
\usage{
read_pcpt_synthetic(path = NULL)

prepare_pcpt_linkage_data(data, change_points = c(2373, 2738))

compare_pcpt_methods(
  data,
  time_dependent = TRUE,
  cox_covariates = NULL,
  linkage_covariates = c(
    "PSA", "racecode", "bmir", "famhx", "txg", "agebase",
    "EDUCATN", "MARISTAT", "EMPLSTAT", "PROFESN"
  )
)
}
\arguments{
\item{path}{Optional path to a compatible space-delimited dataset. By default,
reads the synthetic dataset bundled with the package.}
\item{data}{For \code{prepare_pcpt_linkage_data()}, a subject-level PCPT-format
data frame. For \code{compare_pcpt_methods()}, a list returned by
\code{prepare_pcpt_linkage_data()}.}
\item{change_points}{Two increasing positive follow-up times defining
treatment-effect change points. Defaults are 6.5 and 7.5 years in days.}
\item{time_dependent}{Whether to include the treatment-effect change terms
\code{T1} and \code{T2}.}
\item{cox_covariates}{Optional character vector of Cox-model covariates.}
\item{linkage_covariates}{Character vector of linkage-model covariates.}
}
\value{
\code{read_pcpt_synthetic()} returns a data frame.
\code{prepare_pcpt_linkage_data()} returns a list containing \code{df} and
\code{tdata}. \code{compare_pcpt_methods()} returns estimates, standard errors,
confidence intervals, and fitted objects for IPLW, CC, CC+, and NLAC.
}
\description{
Read the bundled synthetic Prostate Cancer Prevention Trial data, prepare its
counting-process representation, and compare the paper's linkage methods.
The bundled data contain no original participant records.
}
\examples{
pcpt <- read_pcpt_synthetic()
head(pcpt)

\donttest{
prepared <- prepare_pcpt_linkage_data(pcpt)
results <- compare_pcpt_methods(prepared)
results$beta_hat
}
}
