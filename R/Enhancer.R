## This is the proposed enhancing function
GEDE <- function(Y, Est="auto", covariates=NULL, predictors=seq(1:ncol(Y)), HD=FALSE, HD.iter=5, nMAD=2, ...) {
  if (identical(Est,"auto")) Est <- RobEst(Y, covariates=covariates, HD=HD, HD.iter=HD.iter, nMAD=nMAD, ...)
  Tk <- Est$Tk; Lk <- Est$Lk
  sigma2 <- Est$sigma2; K <- Est$K; n <- nrow(Y)
  ## compute the expected values for the test data (Y)
  betahat <- Est$betahat
  if (is.null(covariates)) { #only use the intercept
    mumat <- rep(1,n)%*%t(betahat[1,])
  } else {
    mumat <- cbind(1, covariates)%*%betahat
  }
  ## outliers should be defined by the centered test data
  out.idx <- Hampel(Y-mumat, nMAD=nMAD)
  Y.out <- Y; Y.out[out.idx] <- NA
  suppressWarnings(Y.imputed <- EigenImpute(Est, Y.out, covariates=covariates, predictors=predictors, HD=HD, HD.iter=HD.iter))
  ## Now conduct enhancement based on auto prediction.
  if (sigma2==0) {  #no measurement error
    Ystar <- Y
  } else if (K==0) { #there is no useful, sample-specific signal
    Ystar <- mumat
  } else { #the main case
    if (identical(predictors, seq(1:ncol(Y)))) { #auto-prediction
      Yc <- Y.imputed-mumat
      Tktilde <- sweep(Tk, 2, Lk/(sigma2+Lk), "*")
      Ystar <- mumat +(Yc%*%Tk)%*%t(Tktilde)
    } else { #using only a subset of predictors
      Tk2 <- Tk[predictors,]; Y2 <- Y.imputed[,predictors]
      ss <- svd(sweep(Tk2, 2, sqrt(Lk), "*"))
      U <- ss$u; V <- ss$v; dd <- ss$d
      Tktilde <- sweep(Tk, 2, sqrt(Lk), "*")%*%V
      Utilde <- sweep(U, 2, dd/(sigma2+dd^2), "*")
      Y2c <- Y2-mumat[,predictors]
      Ystar <- mumat +(Y2c%*%Utilde)%*%t(Tktilde)
    }
  }
  dimnames(Ystar) <- dimnames(Y)
  ## Add Ystar and out.idx to Est and return Est
  Est$Ystar <- Ystar; Est$out.idx <- out.idx; Est$NA.idx <- which(is.na(Y))
  return(Est)
}

## A consistent enhancing interface for several methods, using train
## and test data
Enhancer <- function(train, test, covariates.train=NULL, covariates.test=NULL, method=c("GEDE", "lasso", "lasso2"), predictors=seq(1:ncol(train)), mc.cores=2, ...) {
  method <- match.arg(method)
  train <- as.matrix(train); test <- as.matrix(test)
  n <- nrow(train); m <- ncol(train)
  if (method=="GEDE"){
    Est <- RobEst(train, covariates=covariates.train, ...)
    rr <- GEDE(test, Est=Est, covariates=covariates.test, predictors=predictors, ...)
  } else if (method=="lasso") {
    Ystar <- Reduce(cbind, mclapply(1:m, function(i) {
      Xi.idx <- setdiff(predictors, i)
      Xtrain <- cbind(train[,Xi.idx], covariates.train)
      mod.i <- cv.glmnet(x=Xtrain, y=train[,i], type.measure="mse",
                         alpha=1, nfolds=5)
      Xtest <- cbind(test[,Xi.idx], covariates.test)
      predict(mod.i, s=mod.i$lambda.min, newx=Xtest)
    }, mc.cores=mc.cores))
    rr <- list(Ystar=Ystar)
  } else if (method=="lasso2"){ #a "global" CV
    Ystar <- Reduce(cbind, mclapply(1:m, function(i) {
      Xi.idx <- setdiff(predictors, i)
      Xtrain <- cbind(train[,Xi.idx], covariates.train)
      mod.i <- glmnet(x=Xtrain, y=train[,i])
      ## Remarks: (a) the deviance of a Gaussian model is its RSS; (b)
      ## nrow(train)-mod.i$df is (approximately) the residual degrees of
      ## freedom; (c) gcv.i actually equals GCV/n
      gcv.i <- deviance(mod.i)/(n-mod.i$df)^2
      best.lambda <- mod.i$lambda[which.min(gcv.i)]
      Xtest <- cbind(test[,Xi.idx], covariates.test)
      predict(mod.i, s=best.lambda, newx=Xtest)
    }, mc.cores=mc.cores))
    rr <- list(Ystar=Ystar)
  } else {
    stop(paste0("The method you specified, ", method, ", has not been implemented in Enhancer() yet."))
  }
  return(rr)
}
