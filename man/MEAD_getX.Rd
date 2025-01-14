% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{MEAD_getX}
\alias{MEAD_getX}
\title{Prepare reference matrix and its variance for deconvolution}
\usage{
MEAD_getX(ref, cell_types, individuals)
}
\arguments{
\item{ref}{returned by MEAD_preprocessing}

\item{cell_types}{cell types to be considered}

\item{individuals}{subject in reference dataset to use}
}
\value{
a list of:
 *X: reference matrix, G by K
 *V: variance of reference matrix, G by K^2
 *w: a vector of gene weights, length G
}
\description{
Prepare reference matrix and its variance for deconvolution
}
