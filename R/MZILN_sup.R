## Supplemental functions for MZILN

# Overall function -------------------------------------------#####
#' This function performs per sample taxonomic transformation 
#' @description For a sample, we find the log-ratio transformation of all non-zero taxa as well as the design matrix 
#' @param main A sample from the taxonomic abundance matrix
#' @param n.tax Number of total taxa 
MZILN.transformation <- function(df.main, df.covar, verbose = F){
  if (verbose == T){
    print("Begin transformations...") 
  }
  n.tax <- ncol(df.main)
  n.subj <- nrow(df.main)
  outcome <- c()
  predictor <- c()
  for (i in 1:n.subj){
    if (verbose == T){
      if (i %% 50 == 0){
        print(paste("Currently at subject",i,"of",n.subj,"subjects"))
      }
    }
    tax.trans <- tax.transformation(df.main[i,], n.tax)
    log.trans <- tax.trans$log
    matrix.a <- tax.trans$mat
    rownames(matrix.a) <- NULL
    omega <- solve(matrix.a %*% t(matrix.a)) # assuming identity
    matrix.x <- covar.transformation(df.covar[i,], n.tax)
    outcome <- rbind(outcome,expm::sqrtm(omega) %*% log.trans)
    predictor <- rbind(predictor,expm::sqrtm(omega) %*% matrix.a %*% matrix.x)
    
  }
  result <- list(outcome = outcome, predictor = predictor)
  return(result)
}

# Tax transform -------------------------------------------####
#' This function performs per sample taxonomic transformation 
#' @description For a sample, we find the log-ratio transformation of all non-zero taxa as well as the design matrix 
#' @param main A sample from the taxonomic abundance matrix
#' @param n.tax Number of total taxa 
tax.transformation <- function(main, n.tax){
  main <- as.numeric(main)
  nonzero <- which(main != 0) #index of nonzero elements length of this is L. 
  
  # log transformed outcome 
  log.trans <- log(main[nonzero][-length(nonzero)]/main[nonzero][length(nonzero)])
  
  #matrix a
  mat <- matrix(0,ncol = n.tax - 1, nrow = length(nonzero) - 1) # this matrix is of dimensions (L - 1) x K
  for (i in 1:(length(nonzero) - 1)){ # remove the last element
    mat[i,nonzero[i]] <- 1
  }  
  if(nonzero[length(nonzero)] != n.tax){ # check if last non-zero element is last taxa
    mat[,nonzero[length(nonzero)]] <- rep(-1,length(nonzero) - 1)
  } 
  return(list(mat = mat, log = log.trans))
}


# Covar transform -------------------------------------------#####
#' Getting the transformed covariance matrix for one sample 
#' @description The covariance matrix for each sample is the kronecker product between the identity matrix
#'    with K dimensions (number of taxa - 1) and the covariance row vector (added 1 for intercept)
#' @param covar row covariance vector 
#' @param n.tax Number of taxa. This is K + 1
covar.transformation <- function(covar, n.tax){
  covar <- t(matrix(c(1, as.numeric(covar))))
  diag_mat <- diag(x = 1, nrow = n.tax - 1)
  matrix.x <- kronecker(diag_mat, covar)
  return(matrix.x)
}




# Cross validation function for picasso -------------------------------------------####

cvPicasso <- function(x, y, nfolds=10, zeroSDCut=10^(-6), lambda.min.ratio=10^(-3), nLam=100, method="mcp", standardize=F, seed = 11){
  print("Start cross validation...")
  # check if there are x variables with zero variance
  # check small variances of x
  sdX <- apply(x,2,sd)
  xWithNearZeroSd <- which(sdX < zeroSDCut)
  rm(sdX)
  nObsAll <- length(y)
  nearZeroSd <- length(xWithNearZeroSd)
  print(paste("Number of x with near-zero sd: ", length(xWithNearZeroSd)))
  
  # fill those zero-variance x variables with randomly generated standard normal
  if (length(xWithNearZeroSd) > 0){  
    randomNorm <- rnorm(n = length(y)*length(xWithNearZeroSd))
    for (i in 1:length(xWithNearZeroSd)){
      x[,xWithNearZeroSd[i]] <- randomNorm[(1+(i-1)*length(y)):(i*nObsAll)]
    }
    
    # recheck if there are x variables with zero variance after the filling
    sdX.r=apply(x,2,sd)
    xWithNearZeroSd.r=which(sdX.r<zeroSDCut)
    rm(sdX.r)
    
    nearZeroSd.r=length(xWithNearZeroSd.r)
    rm(xWithNearZeroSd.r)
    
    print(paste("Number of x with near-zero sd after correction: ", nearZeroSd.r))
  }

  # calculate lambda max
  lamMax=max(abs(colSums(x*y)))/nObsAll
  print(paste("Lambda max: ",lamMax))
  lamMin=lambda.min.ratio*lamMax
  lamList=seq(lamMax,lamMin,length=nLam)
  print(paste("Lambdas calculated for cross validation:"))
  print(lamList)
  
  # partition the data randomly into nfolds partitions
  parSize=floor(nObsAll/nfolds)
  
  # set seed so that the random partition is reproducible
  set.seed(seed)
  randomShuf=sample(nObsAll, nObsAll)
  sampleInd=list()
  for (i in 1:nfolds){
    if(i<nfolds){
      sampleInd[[i]]=randomShuf[(1+(i-1)*parSize):(i*parSize)]
    }else{
      sampleInd[[i]]=randomShuf[(1+(i-1)*parSize):nObsAll]
    }
  }
  
  # cross validation
  SSE=rep(0,nLam)
  
  print(paste("Initialize sum of square error:"))
  
  #sink(paste("picassoPredict",method,".txt",sep="")) # to avoid output from the predict() function
  
  for(i in 1:nfolds){
    # check if there is zero-variance x in the partitions
    xi=x[-sampleInd[[i]],]
    nObs.i=nrow(xi)
    sdX.i=apply(xi,2,sd)
    xWithNearZeroSd.i=which(sdX.i<zeroSDCut)
    rm(sdX.i)
    
    nearZeroSd.i=length(xWithNearZeroSd.i)
    print(paste("Number of x in partition with near-zero sd: ", nearZeroSd.i))
    
    # fill those zero-variance x in the partition with randomly generated normal
    if (nearZeroSd.i>0){  
      randomNorm=rnorm(n=nObs.i*nearZeroSd.i)
      for (j in 1:nearZeroSd.i){
        xi[,xWithNearZeroSd.i[j]]=randomNorm[(1+(j-1)*nObs.i):(j*nObs.i)]
      }
      
      # recheck if there are x variables with zero variance after the filling
      sdX.i.r=apply(xi,2,sd)
      xWithNearZeroSd.i.r=which(sdX.i.r<zeroSDCut)
      rm(sdX.i.r)
      
      nearZeroSd.i.r=length(xWithNearZeroSd.i.r)
      rm(xWithNearZeroSd.i.r)
      
      print(paste("Number of x in partition with near-zero sd after correction: ", nearZeroSd.i.r))
    }else {
      print("No x with near-zero variance in the partition")
    }
    
    yi=y[-sampleInd[[i]]]
    nObsYi=length(yi)
    nObsXi=nrow(xi)
    print(paste("number of x in partition after zero-variance correction: ",nObsXi))
    print(paste("number of y in partition: ",nObsYi))
    
    colsumYi=sum(abs(yi))/nObs.i
    colsumXiMax=max(abs(colSums(xi)))/nObs.i
    
    print(paste("colsumYi: ",colsumYi))
    print(paste("colsumXiMax: ",colsumXiMax))
    
    print(paste("start cross validation for fold: ",i))
    
    cv.i=picasso::picasso(X=xi,Y=yi,lambda=lamList,method=method,standardize=standardize)
    rm(xi,yi)
    
    if(nearZeroSd.i>0) {
      cv.i$beta[xWithNearZeroSd.i,]=0  
    }
    rm(xWithNearZeroSd.i)
    
    #print("beta in cv.i: ")
    #print(cv.i$beta)
    
    print("Penalized reggression done, start estimating prediction")
    
    testSetSize=length(sampleInd[[i]])
    yHatTest.i=x[sampleInd[[i]],]%*%(cv.i$beta)+matrix(rep(cv.i$intercept,each=testSetSize),nrow=testSetSize)
    rm(cv.i)
    
    resiVecs.i=(yHatTest.i-y[sampleInd[[i]]])
    rm(yHatTest.i)
    SSEi=apply(resiVecs.i,2,function(x)sum(x^2))
    rm(resiVecs.i)
    SSE=SSE+SSEi
    rm(SSEi)
    print(paste("Cross validation done for fold: ",i))
  }
  optiLamInd=which(SSE==min(SSE))
  optiLam=lamList[optiLamInd]
  print(paste("Lambda selected by cross validation",optiLam))
  plot(x = lamList, y = SSE, ylab = "SSE", xlab = "Lambda", pch = 18)
  abline(v = optiLam, col = "red", lwd = 3)
  return(list(optiLam = optiLam, SSE = SSE, lamList = lamList))
}

