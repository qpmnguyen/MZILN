#### Main functions ------------------------------------------------------####

####ZILN Main function ---------------------------------------------------####

#main ZILN function ----------------------------------------------------####
#' @export
#'
ZILN.main <- function(df.main, df.covar, covariates=colnames(df.covar), n.lam = 200, lam.min.ratio = 0.0001,
                      reg.method = 'mcp', n.folds = 5){
  log.data.full <- c() #initializing a main data function
  matrix.data.full <- c() #initializing a main matrix function
  n.sub <- nrow(df.main) #number of subjects



  for (m in (1:n)){
    nonzero <- which(df.main[m,] != 0)
    log.trans <- log.trans.generation(df.main,m,nonzero)
    matrix.a <- matrix.a.generation(nonzero)
    matrix.x <- matrix.x.generation(df.covar,m)
    omega <- solve(matrix.a %*% t(matrix.a))
    log.data.full <- rbind(log.data.full, expm::sqrtm(omega) %*% log.trans)
    matrix.data.full <- rbind(matrix.data.full, expm::sqrtm(omega) %*% matrix.a %*% matrix.x)
  }

  log.data.full <- as.numeric(log.data.full)
  picasso.data <- picasso::picasso(X = matrix.data.full,Y = log.data.full,
                          lambda = NULL, nlambda = n.lam, family = 'gaussian',
                          method = reg.method,type.gaussian = 'naive',
                          standardize = FALSE, verbose = FALSE,lambda.min.ratio = lam.min.ratio)
  lambda.values <- picasso.data$lambda

  ###### cross validation #####
  shuffle <- sample(nrow(matrix.data.full)) #creating a shuffle pattern
  log.data.full <- log.data.full[shuffle] #shuffle this data set
  matrix.data.full <- matrix.data.full[shuffle,] #shuffle this data set similar to the one before
  folds <- cut(seq(1:nrow(matrix.data.full)),breaks=n.folds,labels=FALSE)
  lambda.error <- matrix(nrow = n, ncol = n.folds)

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
    sink(file = NULL)
    predict <- predict.gaussian(picasso.train, matrix.test.full,
                                lambda.idx = c(1:200), y.pred.idx = c(1:200))
    closeAllConnections()
    for (l in (1:n)){
      lambda.error[l,i] <- (sum((predict[,l] - log.test.full)^2))/nrow(log.test.full)
    }
  }

  #lambda errors
  lambda.error <- as.data.frame(lambda.error)
  lambda.error$mean <- rowSums(lambda.error)/5
  lambda.error$lambda <- lambda.values
  v <- which(lambda.error$mean == min(lambda.error$mean)) #lambda value with the lowest

  #beta values
  beta <- as.vector(picasso.data$beta[,v])
  #generating the results

  result <- split(beta, ceiling(seq_along(beta)/))

}


#### end ZILN main function ----------------------------------------------####

#### end main functions --------------------------------------------------####
