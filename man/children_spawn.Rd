% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICRegress.R
\name{children_spawn}
\alias{children_spawn}
\title{Create New Children}
\usage{
children_spawn(
  n_children,
  dna_length = ncol(newdf),
  col_names = c("intercept", colnames(newdf)[2:length(colnames(newdf))]),
  normal = T,
  this_sd = 0.5,
  this_range = c(-10, 10)
)
}
\description{
Generate random population (coefficients) for the genetic algorithm ('DNA' or  'children')
}
