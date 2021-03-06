% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting_functions.R
\name{format_axis_log10}
\alias{format_axis_log10}
\title{Generate axis breaks and labels most suitable 
for log10 transformed variables}
\usage{
format_axis_log10(
  range,
  add_minus_sign = FALSE,
  minor_ticks = 2:9,
  n_labels_skip = 0
)
}
\arguments{
\item{range}{range of values (log10 transformed) for which the breaks 
and labels will be generated}

\item{add_minus_sign}{logical. Add a minus sign to log10 transformed values 
(useful for volcano plots where p-values appear in decreasing order)}

\item{minor_ticks}{set of multiplicative factors used to generate 
minor ticks (set to 2:9 by default). Ignored if NULL.}

\item{n_labels_skip}{Number of labels skipped between two displayed labels.}
}
\value{
a list with elements
\itemize{
\item{breaks}{numeric vector of breaks}
\item{labels}{vector of labels (as expressions)}  
}
}
\description{
Generate axis breaks and labels most suitable 
for log10 transformed variables
}
