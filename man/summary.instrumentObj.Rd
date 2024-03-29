% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.instrumentObj}
\alias{summary.instrumentObj}
\title{Summarize a theta2 object}
\usage{
\method{summary}{instrumentObj}(object, pars = "default", probs = c(0.025, 0.5, 0.975), ...)
}
\arguments{
\item{object}{theta2Obj model object}

\item{pars}{"default" will give the default parameters of the model.
Other options will give only the specified vector of parameter names.
For example, pars = c("theta", "alpha").}

\item{probs}{probability vector for quantiles of posterior estiamtes}

\item{...}{not used}
}
\value{
a \code{data.table} of draws
}
\description{
Produce summary data table for theta2 output.
}
