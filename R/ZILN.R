#### Main functions ------------------------------------------------------####

####ZILN Main function ---------------------------------------------------####

#main ZILN function ----------------------------------------------------####
#' ZILN (Zero inflated Logistic Normal Regression)
#'
#' @export
#'
#' @description ZILN is based upon the Zero Inflated Logistic Normal Regression
#' for microbiome relative abundance data prepared by Dr. Zhigang Li.
#' @param  df.main Data frame of patient relative abundance data
#' @param  df.covar Data frame of patient covariates
#' @param  covariates Specific covariates of interest. The default is all covariates
#'     in \code{df.covar}
#' @param  n.lam Number of \code{lambda} values to control the regularization as documented
#'     in \code{PICASSO}. The default is 200.
#' @param lam.min.ratio Ratio between the upper and lower bound of generated \code{lambda} values.
#'     The default is 0.0001.
#' @param reg.method The method of regularization. The default is \code{mcp}. Other options
#'     include \code{lasso} and \code{scad}
#' @param n.folds The number of cross-validation folds. The default is 5.
#'
#' @return A list containing: \code{results} - a data frame of regression coefficients for each
#'     covariate of interest (as columns) per taxa (as rows); and \code{crossvalidation} - a data
#'     frame of each \code{lambda} value and the mean MSE of each value across cross validation folds

ZILN.main <- function(df.main, df.covar, covariates=colnames(df.covar), n.lam = 200,
                      lam.min.ratio = 0.0001,
                      reg.method = 'mcp', n.folds = 5){

  log.data.full <- c() #initializing a main data function
  matrix.data.full <- c() #initializing a main matrix function
  n.sub <- nrow(df.main) #number of subjects
  #creating a localized version of the covariates data set depending upon which
  #covariates to be selected
  df.covariates <- covar.restriction(df.covar,covariates)

  #transformation of the data
  for (m in (1:n.sub)){
    nonzero <- which(df.main[m,] != 0)
    log.trans <- log.trans.generation(df.main,m,nonzero)
    matrix.a <- matrix.a.generation(nonzero)
    matrix.x <- matrix.x.generation(df.covariates,m)
    omega <- solve(matrix.a %*% t(matrix.a))
    log.data.full <- rbind(log.data.full, expm::sqrtm(omega) %*% log.trans)
    matrix.data.full <- rbind(matrix.data.full, expm::sqrtm(omega) %*% matrix.a %*% matrix.x)
  }
  log.data.full <- as.numeric(log.data.full)

  #picasso for the whole data.
  picasso.data <- picasso::picasso(X = matrix.data.full,Y = log.data.full,
                          lambda = NULL, nlambda = n.lam, family = 'gaussian',
                          method = reg.method,type.gaussian = 'naive',
                          standardize = FALSE, verbose = FALSE,lambda.min.ratio = lam.min.ratio)

  lambda.values <- picasso.data$lambda #extracting lambda values

  ###### cross validation #####
  shuffle <- sample(nrow(matrix.data.full)) #creating a shuffle pattern
  log.data.full <- log.data.full[shuffle] #shuffle the log transformed data set.
  matrix.data.full <- matrix.data.full[shuffle,] #shuffle matrix x similar to the one before
  folds <- cut(seq(1:nrow(matrix.data.full)),breaks=n.folds,labels=FALSE)
  lambda.error <- matrix(nrow = n.lam, ncol = n.folds)
  for(i in (1:n.folds)){
    #Segement your data by fold using the which() function
    testIndexes <- which(folds==i,arr.ind=TRUE)
    test.log <- log.data.full[testIndexes]
    test.covar <- matrix.data.full[testIndexes,]
    train.log <- log.data.full[-testIndexes]
    train.covar <- matrix.data.full[-testIndexes,]
    #set it to  lambda.values
    picasso.train <- picasso::picasso(X = train.covar,Y = train.log,
                             lambda = lambda.values, family = 'gaussian',
                             method = reg.method, type.gaussian = 'naive',
                             standardize = FALSE, verbose = FALSE)
    f <- file()
    sink(file = f)
    predict <- picasso::predict.gaussian(picasso.train, test.covar,
                                lambda.idx = c(1:n.lam), y.pred.idx = c(1:n.lam))
    sink()
    close(f)
    for (l in (1:n.lam)){
      lambda.error[l,i] <- (sum((predict[,l] - test.log)^2))/nrow(as.matrix(test.log))
    }
  }

  #lambda errors
  lambda.error <- as.data.frame(lambda.error)
  lambda.error$mean <- rowSums(lambda.error)/n.folds
  lambda.error$lambda <- lambda.values
  v <- which(lambda.error$mean == min(lambda.error$mean)) #lambda value with the lowest
  #beta values
  beta <- as.vector(picasso.data$beta[,v])
  #generating the results
  list <- list() #empty list
  #result
  result <- split(beta, ceiling(seq_along(beta)/(ncol(df.covariates) + 1)))

  index <- c()
  for (i in (1:length(result))){
    if (isTRUE(all.equal(result[[i]][-1],rep(0,length(covariates))))){
      index[i] <- i
    } else
      index[i] <- NA
  }
  index <- na.omit(index)
  result <- result[-index]
  result <- as.data.frame(do.call("rbind",result))
  rownames(result) <- colnames(df.main)[-c(index,length(colnames(df.main)))]
  colnames(result) <- c("Intercept",covariates)

  list[[1]] <- result
  list[[2]] <- lambda.error[,c("lambda","mean")]
  names(list) <- c("results","crossvalidation")
  return(list)
}


#### end ZILN main function ----------------------------------------------####

#### end main functions --------------------------------------------------####
