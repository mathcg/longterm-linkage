\name{longterm_fit_methods}
\alias{longterm_fit_methods}
\alias{print.longterm_fit}
\alias{coef.longterm_fit}
\alias{vcov.longterm_fit}
\alias{confint.longterm_fit}
\title{Methods for Fitted Long-Term Linkage Models}
\usage{
\method{print}{longterm_fit}(x, ...)

\method{coef}{longterm_fit}(object, ...)

\method{vcov}{longterm_fit}(object, ...)

\method{confint}{longterm_fit}(object, parm, level = 0.95, ...)
}
\arguments{
\item{x, object}{A \code{longterm_fit} object.}
\item{parm}{Coefficients to include, by position or name.}
\item{level}{Confidence level.}
\item{...}{Additional arguments, currently unused.}
}
\value{
\code{coef()} returns a named numeric vector, \code{vcov()} returns a matrix,
and \code{confint()} returns a two-column matrix. \code{print()} returns its
input invisibly.
}
\description{
Extract and display results returned by \code{fit_linkage_cox()}.
}
