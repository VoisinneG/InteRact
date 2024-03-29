% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting_functions.R
\name{plot_volcanos}
\alias{plot_volcanos}
\title{Plot protein enrichement fold-change versus p-value}
\usage{
plot_volcanos(
  res = NULL,
  data = NULL,
  names = NULL,
  idx = NULL,
  labels = NULL,
  show_all_above_thresh = FALSE,
  N_print = 15,
  conditions = NULL,
  p_val_thresh = 0.05,
  fold_change_thresh = 2,
  x0 = NULL,
  c = NULL,
  save_file = NULL,
  xlim = NULL,
  ylim = NULL,
  norm = FALSE,
  both_sides = FALSE,
  show_thresholds = TRUE,
  alpha_segment = 0.2,
  theme_name = "theme_bw",
  size = 1,
  alpha = 1,
  color = rgb(0.75, 0.75, 0.75),
  label_size = 3,
  n_character_max = 8,
  n_labels_skip_x = 0,
  n_labels_skip_y = 1,
  minor_ticks = 2:9,
  ...
)
}
\arguments{
\item{res}{an \code{InteRactome}}

\item{data}{data.frame with columns 'names', 'p_val' and 'fold_change'}

\item{names}{Names of proteins highlighted}

\item{idx}{indices of proteins highlighted (superseeds \code{names})}

\item{labels}{labels for proteins in plot. Must the same length as 
\code{res$names} or \code{data$names}}

\item{show_all_above_thresh}{logical. Highlight all proteins above threshold?}

\item{N_print}{maximum of protein labels to display}

\item{conditions}{conditions to plot}

\item{p_val_thresh}{threshold on p-value to display}

\item{fold_change_thresh}{threshold on fold-change to display}

\item{x0}{parameters x0 of the line dividing the volcano plot according to \code{f(x) = c / (|x|-x0)}.
Ignored unless parameters \code{p_val_thresh} and \code{fold_change_thresh} are set to \code{NULL}}

\item{c}{parameters c of the line dividing the volcano plot according to \code{f(x) = c / (|x|-x0)}.
Ignored unless parameters \code{p_val_thresh} and \code{fold_change_thresh} are set to \code{NULL}}

\item{save_file}{path of output file (.pdf)}

\item{xlim}{range of x values}

\item{ylim}{range of y values}

\item{norm}{Use normalized fold-changes}

\item{both_sides}{logical. Shading on right and left upper graphs.}

\item{show_thresholds}{Show thresholds using red lines?}

\item{alpha_segment}{transparency of threshold segments}

\item{theme_name}{name of the ggplot2 theme function to use ('theme_gray' by default)}

\item{size}{dot size}

\item{alpha}{dot transparency}

\item{color}{dot color}

\item{label_size}{size of labels (5 by default)}

\item{n_character_max}{max number of label characters}

\item{n_labels_skip_x}{Number of labels skipped between two displayed labels on x axis.}

\item{n_labels_skip_y}{Number of labels skipped between two displayed labels on y axis.}

\item{minor_ticks}{set of multiplicative factors used to generate 
minor ticks (set to 2:9 by default). Ignored if NULL.}

\item{...}{parameters passed to \code{geom_text_repel()}}
}
\value{
a plot
}
\description{
Plot protein enrichement fold-change versus p-value
}
