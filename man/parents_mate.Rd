% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICRegress.R
\name{parents_mate}
\alias{parents_mate}
\title{Genetic Crossover}
\usage{
parents_mate(
  dna_set = dna_pool,
  first_parent_indices,
  second_parent_indices,
  mutation_rate,
  random_mate = T,
  mutation_granularity = 1,
  percent_explained = solved,
  mode = "continuous"
)
}
\description{
Creates new children based on the current successful dnapool
}
