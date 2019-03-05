## Supplemental functions for MZILN

tax.transformation <- function(main, n.tax, type){
  main <- as.numeric(main)
  nonzero <- which(main != 0) #index of nonzero elements
  if(type == "log"){
    log.trans <- log(main[nonzero][-length(nonzero)]/main[nonzero][length(nonzero)])
    return(log.trains)
  }
  if (type == "matrix.a"){
    row <- rep(0, n.tax - 1)
    if(nonzero[length(nonzero)] == n.tax){
      row[nonzero[-length(nonzero)]] <- 1
    } else {
      row[nonzero[-length(nonzero)]] <- 1
      row[nonzero[length(nonzero)]] <- -1
    }
    return(row)
  }
}

covar.transformation <- function(covar, n.tax){
  covar <- c(1, as.numeric(covar))
  list <- rep(list(t(covar)), n.tax - 1)
  matrix.x <- as.matrix(Matrix::bdiag(list))
  return(matrix.x)
}




