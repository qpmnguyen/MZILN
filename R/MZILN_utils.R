###-utilities functions that are used by the package but not exported ----------####


nsubj <- function(data.main){
  return(nrow(data.main))
}

idsubj <- function(data.main){
  return(rownames(data.main))
}

covar.name <- function(data.covar){
  return(colnames(data.covar))
}

ntaxa <- function(data.main){
  return(ncol(data.main))
}


taxa.name <- function(data.main){
  return(colnames(data.main))
}

# check for all zero columns and all zero rows
zero.check <- function(data.main){
  check_rsum <- rowSums(data.main)
  check_csum <- colSums(data.main)
  rzero <- which(check_rsum == 0)
  czero <- which(check_csum == 0)
  return(list(rzero = rzero, czero = czero))
}
### end of utilities functions -------------------------------------------------####
