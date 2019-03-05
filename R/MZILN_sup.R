## Supplemental functions for MZILN

MZILN.transformation <- function(df.main, df.covar){
  n.tax <- ncol(df.main)
  n.subj <- nrow(df.main)
  
  # getting log transformations  
  log.trans <- apply(df.main,1, tax.transformation, n.tax, "log")
  log.trans <- unlist(log.trans, use.names = F) # convert list into one big vector
  
  # getting matrix.a
  matrix.a <- apply(df.main,1, tax.transformation, n.tax, "matrix.a")
  matrix.a <- t(matrix.a)
  rownames(matrix.a) <- NULL # removing uncessary rownames
  
  # getting matrix.x
  matrix.x <- c()
  for (i in 1:nsubj){
    matrix.x <- rbind(matrix.x,covar.transformation(df.covar[i,], n.tax))
  }
  result <- list(log.trans = log.trans, matrix.a = matrix.a, matrix.x = matrix.x)
}


n.tax <- ncol(test.main)
test <- apply(test.main, 1, tax.transformation, n.tax, "log")
test.2 <- apply(test.main, 1, tax.transformation, n.tax, "matrix.a")
test.3 <- c()
for (i in 1:nrow(test.main)){
  test.3 <- cbind(test.3,covar.transformation(test.covariates[i,],n.tax))
}

test.3 <- lapply(test.covariates,covar.transformation, n.tax)
length(test.3)

test.3[1:10,1:10]


tax.transformation <- function(main, n.tax, type){
  main <- as.numeric(main)
  nonzero <- which(main != 0) #index of nonzero elements
  if(type == "log"){
    log.trans <- log(main[nonzero][-length(nonzero)]/main[nonzero][length(nonzero)])
    return(log.trans)
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

test.4 <- covar.transformation(test.covariates[1,], n.tax)
dim(test.4)
test.4[1:10,1:10]
