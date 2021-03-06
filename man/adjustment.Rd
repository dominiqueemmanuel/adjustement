% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjustment.r
\name{adjustment}
\alias{adjustment}
\title{A function to calculate satistical adjustment weight}
\usage{
adjustment(data, weight = NULL, margins = list(), weight_min = 0.1,
  weight_max = 10, ...)
}
\arguments{
\item{data}{a data frame.}

\item{weight}{a nonnegative numeric vector (with positive sum) of length \code{nrow(data)}.}

\item{margins}{a list specifing the margins constrains.}

\item{weight_min}{a numeric constant greater or equal to 1.}

\item{weight_max}{a nonnegative numeric constant.}
}
\description{
This function allows to calculate satistical adjustment weight given some margins (equalities or inequalities).
}
\examples{
data <- iris
summary(data)
margins <- list(list(var_name = 'Species',value=c(setosa=0.315),
min=c(versicolor=0.45),
max=c(versicolor=0.01)))
w <- adjustment(data = data,margins = margins)$w
data$w<-w
library(dplyr)
data\%>\%
  group_by(Species)\%>\%
  summarise(eff=sum(w))\%>\%
  merge(data\%>\%summarise(eff0=sum(w)))\%>\%
  mutate(freq = eff/eff0)
}
\keyword{adjustment}
\keyword{calibration}
\keyword{weighting}

