% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InteRact.R
\name{append_annotations}
\alias{append_annotations}
\title{Append annotations to an \code{InteRactome}}
\usage{
append_annotations(res, annotations = NULL, name_id = "Protein.IDs",
  organism = "mouse")
}
\arguments{
\item{res}{an \code{InteRactome}}

\item{annotations}{type of annotations to append}

\item{name_id}{column name used to map protein identifiers}

\item{organism}{organism for which the annotations have to be appended}
}
\value{
an \code{InteRactome}
}
\description{
Append annotations to an \code{InteRactome}
}