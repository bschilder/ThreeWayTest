% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/T_eta.R
\name{T_eta}
\alias{T_eta}
\title{T_eta}
\usage{
T_eta(z_vector, cov_mat, eta)
}
\arguments{
\item{z_vector}{Column vectorized data matrix with rows represent 
phenotype and columns represent genotype.}

\item{cov_mat}{Estimated covariance matrix of z_vector.}

\item{eta}{Truncated value.}
}
\value{
A numeric value represent calculated statistic T_eta.
}
\description{
Calculation of T_eta with eta fixed.
}
\examples{
z_vector<-MASS::mvrnorm(1,mu=rep(0,9),Sigma = diag(nrow = 9, ncol = 9))
ThreeWayTest::T_eta(z_vector=z_vector, 
cov_mat=diag(nrow = 9, ncol = 9), eta=0.5)
}
