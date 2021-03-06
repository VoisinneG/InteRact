% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main_functions.R
\name{merge_proteome}
\alias{merge_proteome}
\title{Add protein abundance to an \code{InteRactome}}
\usage{
merge_proteome(res, ...)
}
\arguments{
\item{res}{an \code{InteRactome}}

\item{...}{Parameters passed to \code{proteinRuler::map_proteome()}}
}
\description{
Add protein abundance to an \code{InteRactome}. For multiple identifiers, 
the abundance of the first match in the proteome dataset is returned.
}
