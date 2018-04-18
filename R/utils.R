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



#checking for all zero subjects
zero.sub <- function(data.main){
  all.zero <- c()
  for (i in (1:nrow(data.main))){
    if (all(data.main[i,] == 0)){
      all.zero[i] <- i
    } else {
      all.zero[i] <- NA
    }
  }
  all.zero <- as.vector(stats::na.omit(all.zero))
  return(all.zero)
}

# Checking for all zero taxas
zero.taxa <- function(data.main){
  all.zero <- c()
  for (i in (1:ncol(data.main))){
    if (all(data.main[,i] == 0)){
      all.zero[i] <- i
    } else {
      all.zero[i] <- NA
    }
  }
  all.zero <- as.vector(stats::na.omit(all.zero))
  return(all.zero)
}


# Utilities for ZILN main
# taxa.num is defined locally in the MZILN.main function.
# The functions below should inherit taxa.num variable defined locally.
# this function is used when the last taxon is zero
# The inheritance method is to define an environment which is the parent into which the expression is called.
# this means that the 'local'environment is treated as the default environment for these functions, which
# are defined outside of the parent function but called inside the parernt function.
defining.row.1 <- function(x,y,env = parent.frame()){
  row <- rep(0, (env$taxa.num)-1)
  row[x] <- 1
  row[y] <- -1
  return(row)
}

#when the last taxon is non-zero
defining.row.2 <- function(x, env = parent.frame()){
  row <- rep(0, (env$taxa.num)-1)
  row[x] <- 1
  return(row)
}

# generating the matrix a for one patient
matrix.a.generation <- function(nonzero,env=parent.frame()){
  taxa.num = env$taxa.num
  if (max(nonzero) == env$taxa.num){
    nonzero.loop <- nonzero[-length(nonzero)] #took off the last element
    matrix.a <- matrix(nrow = length(nonzero.loop), ncol = (env$taxa.num)-1)
    for (i in (1:length(nonzero.loop))){
      matrix.a[i,] <- defining.row.2(nonzero.loop[i])
    }

  } else {
    nonzero.loop <- nonzero
    matrix.a <- matrix(nrow = length(nonzero.loop)-1, ncol = (env$taxa.num)-1)
    for (k in (1:(length(nonzero.loop)-1))){
      matrix.a[k,] <- defining.row.1(nonzero.loop[k],nonzero.loop[length(nonzero.loop)])
    }
  }
  return(matrix.a)
}

#generating log transformation for one patient

log.trans.generation <- function(df,j,nonzero){
  nonzero.values <- df[j,nonzero]
  log.trans <- c()
  for (i in (1:(length(nonzero.values)-1))){
    log.trans[i] <- log(nonzero.values[i]/nonzero.values[length(nonzero.values)])
  }
  log.trans <- as.matrix(as.numeric(log.trans))
  return(log.trans)
}

#generating matrix x for one patient
matrix.x.generation <- function(df.covar,j,env=parent.frame()){
  covar <- as.matrix(c(1,as.numeric(as.vector(df.covar[j,]))))
  covar <- t(covar)
  covar.list <- list() #empty list with all the covariance
  covar.list <- rep(list(covar), (env$taxa.num)-1)
  matrix.x <- as.matrix(Matrix::bdiag(covar.list)) #the gigantic diagonal matrix
  return(matrix.x)
}

#Reduce the size of the covariates according to what the user wants.
covar.restriction <- function(df.covar, covariates){
  covar.extract <- df.covar[,covariates]
  return(covar.extract)
}




### end of utilities functions -------------------------------------------------####
