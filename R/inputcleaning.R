


# main functions -------------------------------------------------------------------------------------

#### confirmation of features -----------------------------------------------------####

#' Confirmation of data features
#'
#' @export
#' @param data.main Data frame containing Relative Abundances for each sample.
#' @param data.covar Data frame containing covariance data for each sample.
#' @param covariates Names of covariates of interest.
#'
#' @details ZILN generally takes three main arguments: First, a relative abundance
#'     of all taxa per subject, which is \code{data.main}. Second, a separate
#'     data frame containing all covariates by subject, which is \code{data.covar}.
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
#'     #loading data frames
#'     data(test.main)
#'     data(test.covariates)
#'     features.confirm(test.main, test.covariates)

features.confirm <- function(data.main, data.covar, covariates = covar.name(data.covar)){
  if (is.data.frame(data.main) == FALSE & is.data.frame(data.covar) == FALSE &
      is.vector(covariates) == FALSE){
    return("Your data must be in data frame format and your covariates must be in vector format! \n\n")
  } else {
    cat ("Your data format is correct! \n\n")
    cat(paste("The number of subjects is "), nsubj(data.main), "\n\n")
    if (nsubj(data.main) != nsubj(data.covar)) {
      cat("The number of subjects is inconsistent between data frames \n\n")
    } else {
      cat("The number of subjects is consistent between data frames \n\n")
    }
    cat("First 10 Subject IDs", utils::head(idsubj(data.main), n = 10L), "\n\n")
    cat("Possible Covariates \n", covar.name(data.covar), "\n\n")
    cat("Interested Covariates \n", covariates, "\n\n")
    cat("Number of taxa \n", ntaxa(data.main), "\n\n")
    cat("First 10 Taxa Names \n", utils::head(taxa.name(data.main), n = 10L), "\n\n")
    cat("Number of all zero subjects \n", length(zero.sub(data.main)), "\n\n")
    cat("Number of all zero taxa \n", length(zero.taxa(data.main)))
  }
}
#### end confirmation of features --------------------------------------####

#### Input cleaning ----------------------------------------------------####

#' Input cleaning
#'
#' @export
#' @usage remove.zeros(data.main)
#' @param data.main A data frame containing relative abundance data for all patients.
#' @description  \code{remove.zeros} takes in a microbiome relative abundance data frame and remove
#'    all samples with zeros across all taxa and remove taxa with all zeros across samples
#' @return A data frame with all-zero columns and all-zero rows removed.
#' @examples
#' data(test.main) #loading the sample data included in ZILN
#' remove.zeros(test.main)

remove.zeros <- function(data.main){
  data.main <- data.main[-zero.sub(data.main),-zero.taxa(data.main)]
  return(data.main)
}


#### end input cleaning ------------------------------------------------####




# end of main functions --------------------------------------------------------------
