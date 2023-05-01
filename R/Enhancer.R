## This is the proposed enhancing function
GEDE <- function(Y, Est, predictors=seq(1:ncol(Y)), HD=FALSE, HD.iter=5, nMAD=3, verbose=FALSE, ...) {
  mumat <- Est$mumat; Tk <- Est$Tk; Lk <- Est$Lk
  sigma2 <- Est$sigma2; K <- Est$K; n <- nrow(Y)
  ## outliers should be defined by Y, not the training data
  out.idx <- Hampel(Y, nMAD=nMAD)
  Y.out <- Y; Y.out[out.idx] <- NA
  Y.imputed <- EigenImpute(Est, Y.out, predictors=predictors, HD=HD, HD.iter=HD.iter)
  ## Now conduct enhancement based on auto prediction.
  if (sigma2==0) {  #no measurement error
    Xhat <- Y
  } else if (K==0) { #there is no useful, sample-specific signal
    Xhat <- mumat
  } else { #the main case
    if (identical(predictors, seq(1:ncol(Y)))) { #auto-prediction
      Yc <- Y.imputed-mumat
      Tktilde <- sweep(Tk, 2, Lk/(sigma2+Lk), "*")
      Xhat <- mumat +(Yc%*%Tk)%*%t(Tktilde)
    } else { #using only a subset of predictors
      Tk2 <- Tk[predictors,]; Y2 <- Y.imputed[,predictors]
      ss <- svd(sweep(Tk2, 2, sqrt(Lk), "*"))
      U <- ss$u; V <- ss$v; dd <- ss$d
      Tktilde <- sweep(Tk, 2, sqrt(Lk), "*")%*%V
      Utilde <- sweep(U, 2, dd/(sigma2+dd^2), "*")
      Y2c <- Y2-mumat[,predictors]
      Xhat <- mumat +(Y2c%*%Utilde)%*%t(Tktilde)
    }
  }
  dimnames(Xhat) <- dimnames(Y)
  if (verbose) {
    ## Add Xhat and out.idx to Est and return Est
    Est$Xhat <- Xhat; Est$out.idx <- out.idx; Est$NA.idx <- which(is.na(Y))
    return(Est)
  } else {
    return(Xhat)
  }
}

## A consistent enhancing interface for several methods
Enhancer <- function(train, test, covariates.train=NULL, covariates.test=NULL, method=c("GEDE", "lasso", "lasso2"), predictors=seq(1:ncol(train)), mc.cores=2, ...) {
  method <- match.arg(method)
  train <- as.matrix(train); test <- as.matrix(test)
  n <- nrow(train); p <- ncol(train)
  if (method=="GEDE"){
    Est <- RobEst(train, covariates=covariates.train, ...)
    Xhat <- GEDE(test, Est=Est, predictors=predictors, ...)
  } else if (method=="lasso") {
    Xhat <- Reduce(cbind, mclapply(1:p, function(i) {
      Xi.idx <- setdiff(predictors, i)
      Xtrain <- cbind(train[,Xi.idx], covariates.train)
      mod.i <- cv.glmnet(x=Xtrain, y=train[,i], type.measure="mse",
                         alpha=1, nfolds=5)
      Xtest <- cbind(test[,Xi.idx], covariates.test)
      predict(mod.i, s=mod.i$lambda.min, newx=Xtest)
    }, mc.cores=mc.cores))
  } else if (method=="lasso2"){ #a "global" CV
    Xhat <- Reduce(cbind, mclapply(1:p, function(i) {
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
  } else {
    stop(paste0("The method you specified, ", method, ", has not been implemented in Enhancer() yet."))
  }
  return(Xhat)
}
