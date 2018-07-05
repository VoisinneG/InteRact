% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InteRact.R
\name{order_interactome}
\alias{order_interactome}
\title{Order proteins within an \code{InteRactome}}
\usage{
order_interactome(res, var_p_val = "min_p_val", p_val_breaks = c(1, 0.1,
  0.05, 0.01))
}
\arguments{
\item{res}{an \code{InteRactome}}

\item{var_p_val}{name of the p-value variable}

\item{p_val_breaks}{numeric vector to discretize p-value}
}
\value{
an \code{InteRactome}
}
\description{
Order proteins within an \code{InteRactome}
}