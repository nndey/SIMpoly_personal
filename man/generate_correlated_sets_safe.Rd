% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_marker_intersection.R
\name{generate_correlated_sets_safe}
\alias{generate_correlated_sets_safe}
\title{Safe Correlated Set Membership Data}
\usage{
generate_correlated_sets_safe(
  n = 1000,
  x = 5,
  p = rep(0.3, x),
  rho = 0.5,
  multi_prob = c(0.5, 0.3, 0.2)
)
}
\arguments{
\item{n}{Integer. Total number of elements (e.g., SNPs).}

\item{x}{Integer. Number of sets.}

\item{p}{Numeric vector of length \code{x}. Marginal probability of inclusion.}

\item{rho}{Numeric. Correlation between sets (overlap).}

\item{multi_prob}{Probabilities for assigning "None" elements to 1–3 sets.}
}
\value{
A list: counts of set overlaps and logical membership matrix.
}
\description{
Safer version of \code{generate_correlated_sets()} that avoids sampling errors when x = 1.
}
