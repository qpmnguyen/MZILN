#' Simulated data frame of microbiome relative abundance.
#'
#' A dataset containing the relative abundance of 260 taxa (as columns) for
#' 200 patients. This dataset is simulated by Dr. Zhigang Li.
#'
#' @format A data frame with 200 rows and 260 variables:
#' \describe{
#'   \item{taxa1}{Relative abundance}
#'   \item{taxa2}{Relative abundance}
#'   ...
#' }
"test.main"

#' Simulated data frame of patient covariates accompanying \code{test.main}.
#'
#' A dataset containing covariate data (4 covariates) for
#' 200 patients. This dataset is simulated by Dr. Zhigang Li.
#'
#' @format A data frame with 200 rows and 4 variables:
#' \describe{
#'   \item{covar1}{Continuous}
#'   \item{covar2}{Continuous}
#'   \item{covar3}{Continuous}
#'   \item{covar4}{Discrete}
#' }
"test.covariates"

#' True covariate coefficients accompanying test data set \code{test.main}.
#'
#' A data set containing true coefficient values for all covariates (4 covariates) on
#' taxa (259 taxa). This dataset is simulated by Dr. Zhigang Li.
#'
#' @format A data frame with 259 rows and 5 variables:
#' \describe{
#'   \item{Intercept}{Continuous}
#'   \item{covar1}{Continuous}
#'   \item{covar2}{Continuous}
#'   \item{covar3}{Continuous}
#'   \item{covar4}{Continuous}
#' }
"key"
