


# main functions -------------------------------------------------------------------------------------

#### confirmation of features -----------------------------------------------------####

#' Confirmation of data features
#'
#' @export
#' @param df.main Data frame containing Relative Abundances for each sample.
#' @param df.covar Data frame containing covariance data for each sample.
#' @param covariates Names of covariates of interest. The default is all covariates
#'                   in \code{df.covar}
#'
#' @details MZILN generally takes three main arguments: First, a relative abundance
#'     of all taxa per subject, which is \code{df.main}. Second, a separate
#'     data frame containing all covariates by subject, which is \code{df.covar}.
#'     Third, an optional vector containing the names of a subset of covariates that is
#'     of specific interest \code{covariates}. The default mode for \code{covariates}
#'     is the entire columns of covariates in the data frame. ZILN expects subject ids to
#'     be the rownames of the data frame and each column should represent either a taxa's
#'     RA or a covariate value. This function checks all of those required features, display
#'     them as well as check whether or not the formatting is correct.
#'
#' @return Relevant information about the input printed to the console which include: data format,
#'     concordance between data frames, number of subjects, 10 subject IDs, total input
#'     covariates, covariates of interest, Number of taxa, first 10 taxa names, number of zero subjects,
#'     number of zero taxa.
#' @examples
#' #loading data frames
#' data(test.main)
#' data(test.covariates)
#' features.confirm(test.main, test.covariates, covariates = c("covar1","covar4"))

features.confirm <- function(df.main, df.covar, covariates = covar.name(df.covar)){
  if (is.data.frame(df.main) && is.data.frame(df.covar) &&
      is.vector(covariates)){
    cat("Your data format is correct! \n\n")
  } else {
    stop("Your data must be in data frame format and your covariates must be in vector format!")
  }
  message(paste("The number of subjects is "), nsubj(df.main), "\n\n")
  if (nsubj(df.main) != nsubj(df.covar)) {
    warning("The number of subjects is inconsistent between data frames", call.=TRUE)
  } else {
    message("The number of subjects is consistent between data frames \n\n")
  }
  cat("First 10 Subject IDs", utils::head(idsubj(df.main), n = 10L), "\n\n")
  if (isTRUE(all.equal(covariates %in% covar.name(df.covar), rep(TRUE, length(covariates))))){
    cat("Possible Covariates \n", covar.name(df.covar), "\n\n")
    cat("Interested Covariates \n", covariates, "\n\n")
  } else {
    stop("Your specific covariates does not fit covariates listed in your covariate data frame!")
  }
  cat("Number of taxa \n", ntaxa(df.main), "\n\n")
  cat("First 10 Taxa Names \n", utils::head(taxa.name(df.main), n = 10L), "\n\n")
  cat("Number of all zero subjects \n\n", length(zero.sub(df.main)), "\n\n")
  cat("Number of all zero taxa \n\n", length(zero.taxa(df.main)))
  if (length(zero.sub(df.main)) != 0 | length(zero.taxa(df.main)) != 0){
    message("You should use the remove.zeros function to take care of all zero
              taxa and all zero subjects")
  }
}

#### end confirmation of features --------------------------------------####

#### Input cleaning ----------------------------------------------------####

#' Input cleaning
#'
#' @export
#' @usage remove.zeros(df.main)
#' @param df.main A data frame containing relative abundance data for all patients.
#' @description  \code{remove.zeros} takes in a microbiome relative abundance data frame and remove
#'    all samples with zeros across all taxa and remove taxa with all zeros across samples
#' @return A data frame with all-zero columns and all-zero rows removed.
#' @examples
#' data(test.main) #loading the sample data included in MZILN
#' remove.zeros(test.main)

remove.zeros <- function(df.main){
  df.main <- df.main[-zero.sub(df.main),-zero.taxa(df.main)]
  return(df.main)
}


#### end input cleaning ------------------------------------------------####




# end of main functions --------------------------------------------------------------
