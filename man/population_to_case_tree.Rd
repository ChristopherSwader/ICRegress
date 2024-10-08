% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICRegress.R
\name{population_to_case_tree}
\alias{population_to_case_tree}
\title{Population tree to case tree}
\usage{
population_to_case_tree(
  cutlist,
  population_data,
  training_set,
  testing_set,
  training_indices = NULL
)
}
\arguments{
\item{cutlist}{This is the blueprint for building the agglomeration tree, exported from the agglom_tree() function.}

\item{population_data}{Population-based decision tree data from agglom_tree function.}
}
\description{
Applies population differences (from the agglomerative tree) as case-based differences to data set
}
