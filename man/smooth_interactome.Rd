% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main_functions.R
\name{smooth_interactome}
\alias{smooth_interactome}
\title{Smooth, using a moving average across conditions, selected variables of an \code{InteRactome}}
\usage{
smooth_interactome(
  res,
  n = 1,
  order_conditions = NULL,
  var_smooth = c("fold_change", "p_val")
)
}
\arguments{
\item{res}{an \code{InteRactome}}

\item{n}{integer, radius of the moving avergae (number of points extending on each side of the 
center point on which the average is computed)}

\item{order_conditions}{a numeric vector ordering conditions in \code{res$conditions}}

\item{var_smooth}{variables on which the moving average will be computed}
}
\value{
an smoothed \code{InteRactome}
}
\description{
Smooth, using a moving average across conditions, selected variables of an \code{InteRactome}
}
