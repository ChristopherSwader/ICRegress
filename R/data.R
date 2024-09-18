#' A comparison of ICR to other methods
#'
#' The parameters and results data that emerge from the package's benchmarking()
#' function. Parameters used can be found in the package vignette.
#'
#' @format ## `benchmark_results`
#' A data frame with 720 rows 16 columns:
#' \describe{
#'   \item{method}{Method applied to the same underlying data set and ground assumptions.}
#'   \item{weighted r2}{A weighted R-squared value, where the weights depend on the proportion of overall cases in the separate model subgroups.}
#'   \item{discrete solution}{The discrete solution based on the target solution threshold. An alternative parameter used to initially tune the ICR genetic algorithm.}
#'   \item{case assignment score}{A higher score represents that the found subgroups accurately reflect the underlying ground truth of the fake data.}
#'   \item{model similarity error}{A lower model error indicates that the identified models more closely reflect the underlying ground truth of the fake data's models.}
#'   \item{seed}{A random seed.}
#'   \item{generations}{The number of generations used for training for the ICR method's genetic algorithm. Not relevant for the comparison methods.}
#'   \item{solution_thresh}{The threshold on each side of the solution that is considered discretely solved during the training phase.}
#'   \item{n_dna_strands}{The number of models preserved after each generation of training.}
#'   \item{real_model_n}{Underlying number of true models and subgroups of the underlying ground truth fake data.}
#'   \item{N_IVs}{Underlying number of true model-relevant independent variables of the underlying ground truth fake data.}
#'   \item{error_sd}{Underlying error of the ground truth separate subgroups and models. More error indicates messier data.}
#'   \item{test_set_model_closeness}{Deprecated. This is no longer in use as a parameter.}
#'   \item{n models}{The number of models found by ICR and used for comparability for also cluster regression and mixture regression.}
#'   \item{notes}{A description for why NAs are found in certain results.}
#'   \item{parameter_ID}{A unique ID describing a parameter combination.}
#'   ...
#' }
#' @source ICRegress R Package, benchmarking() function.
"benchmark_results"

#' Data from the 2016 European Social Survey
#'
#' An extract of cross-sectional data from the 2016 European Social Survey data is used in order to demonstrate the use of the ICR Method for educational purposes.
#' Full documentation of the ESS and its data are available on its website: https://ess-search.nsd.no/en/study/f8e11f55-0c14-4ab3-abde-96d3f14d3c76.
#' Included variables are a limited selection of those with a theoretical and empirical relationship to happiness.
#'
#' @format ## `ess2016`
#' A data frame with 44387 rows 8 columns:
#' \describe{
#'   \item{happy}{A higher value means greater happiness: "Taking all things together, how happy would you say you are?"}
#'   \item{female}{Dichotomous, a higher value indicates female.}
#'   \item{age}{Age in years.}
#'   \item{cntry}{Country.}
#'   \item{ppltrst}{a higher value indicates higher generalized trust: "enerally speaking, would you say that most people can be trusted, or that you can’t be too careful in dealing with people? }
#'   \item{sclmeet}{A higher value indicates higher frequency of social meetings: "how often do you meet socially with friends, relatives or work colleagues?" (Never, Less than once a month, Once a month, Several times a month, Once a week, Several times a week, Every day)}
#'   \item{income}{Income in deciles.}
#'   \item{health}{Self-rated health: "How is your health37 in general?"}
#'   ...
#' }
#' @source European Social Survey European Research Infrastructure (ESS ERIC). (2020). ESS8 - integrated file, edition 2.2 [Data set]. Sikt - Norwegian Agency for Shared Services in Education and Research. https://doi.org/10.21338/ESS8E02_2
"ess2016"
