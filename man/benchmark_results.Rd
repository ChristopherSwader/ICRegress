% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{benchmark_results}
\alias{benchmark_results}
\title{A comparison of ICR to other methods}
\format{
\subsection{\code{benchmark_results}}{

A data frame with 720 rows 16 columns:
\describe{
\item{method}{Method applied to the same underlying data set and ground assumptions.}
\item{weighted r2}{A weighted R-squared value, where the weights depend on the proportion of overall cases in the separate model subgroups.}
\item{discrete solution}{The discrete solution based on the target solution threshold. An alternative parameter used to initially tune the ICR genetic algorithm.}
\item{case assignment score}{A higher score represents that the found subgroups accurately reflect the underlying ground truth of the fake data.}
\item{model similarity error}{A lower model error indicates that the identified models more closely reflect the underlying ground truth of the fake data's models.}
\item{seed}{A random seed.}
\item{generations}{The number of generations used for training for the ICR method's genetic algorithm. Not relevant for the comparison methods.}
\item{solution_thresh}{The threshold on each side of the solution that is considered discretely solved during the training phase.}
\item{n_dna_strands}{The number of models preserved after each generation of training.}
\item{real_model_n}{Underlying number of true models and subgroups of the underlying ground truth fake data.}
\item{N_IVs}{Underlying number of true model-relevant independent variables of the underlying ground truth fake data.}
\item{error_sd}{Underlying error of the ground truth separate subgroups and models. More error indicates messier data.}
\item{test_set_model_closeness}{Deprecated. This is no longer in use as a parameter.}
\item{n models}{The number of models found by ICR and used for comparability for also cluster regression and mixture regression.}
\item{notes}{A description for why NAs are found in certain results.}
\item{parameter_ID}{A unique ID describing a parameter combination.}
...
}
}
}
\source{
ICRegress R Package, benchmarking() function.
}
\usage{
benchmark_results
}
\description{
The parameters and results data that emerge from the package's benchmarking()
function. Parameters used can be found in the package vignette.
}
\keyword{datasets}
