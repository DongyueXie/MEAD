% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_R01.R
\name{get_R01}
\alias{get_R01}
\title{Estimate the 0-1 Correlation matrix by multiple testing}
\usage{
get_R01(X, alpha = 0.1)
}
\arguments{
\item{X}{gene by individual matrix for estimating correlations}

\item{alpha}{FDR level for claiming correlations}
}
\value{
a sparse matrix, 1 indicates correlated.
}
\description{
Estimate the 0-1 Correlation matrix by multiple testing
}
\details{
This function implements the one sample Large-Scale Correlation Tests With Normal Approximation
proposed by Cai and Liu(2016)
}
