## This is the proposed enhancing function
GEDE <- function(Y, Est="auto", predictors=seq(1:ncol(Y)), HD=FALSE, HD.iter=5, nMAD=3, verbose=FALSE, ...) {
  if (identical(Est,"auto")) Est <- RobEst(Y, HD=HD, HD.iter=HD.iter, ...)
  muhat <- Est$muhat; Tk <- Est$Tk; Lk <- Est$Lk
  sigma2 <- Est$sigma2; K <- Est$K; n <- nrow(Y)
  ## outliers should be defined by Y, not the training data
  out.idx <- Hampel(Y, nMAD=nMAD)
  Y.out <- Y; Y.out[out.idx] <- NA
  Y.imputed <- EigenImpute(Est, Y.out, HD=HD, HD.iter=HD.iter)
  ## Now conduct enhancement based on auto prediction.
  if (sigma2==0) {  #no measurement error
    Xhat <- Y
  } else if (K==0) { #no useful, sample-specific signal
    Xhat <- rep(1,n)%*%t(muhat)
  } else { #the main case
    if (identical(predictors, seq(1:ncol(Y)))) { #auto-prediction
      ## Ytilde <- sweep(Y.imputed, 2, muhat)%*%Tk
      ## Ltilde <- Lk/(sigma2+Lk)
      ## Xhat <- rep(1,n)%*%t(muhat) +Ytilde%*%t(sweep(Tk, 2, Ltilde, "*"))
      Yc <- sweep(Y.imputed, 2, muhat)
      Tktilde <- sweep(Tk, 2, Lk/(sigma2+Lk), "*")
      Xhat <- rep(1,n)%*%t(muhat) +(Yc%*%Tk)%*%t(Tktilde)
    } else { #using only a subset of predictors
      Tk2 <- Tk[predictors,]; Y2 <- Y.imputed[,predictors]
      ss <- svd(sweep(Tk2, 2, sqrt(Lk), "*"))
      U <- ss$u; V <- ss$v; dd <- ss$d
      Tktilde <- sweep(Tk, 2, sqrt(Lk), "*")%*%V
      Utilde <- sweep(U, 2, dd/(sigma2+dd^2), "*")
      Y2c <- sweep(Y2, 2, muhat[predictors])
      Xhat <- rep(1,n)%*%t(muhat) +(Y2c%*%Utilde)%*%t(Tktilde)
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
Enhancer <- function(train, test, method=c("GEDE", "lasso"), predictors=seq(1:ncol(train)), ...) {
  method <- match.arg(method)
  train <- as.matrix(train); test <- as.matrix(test)
  n <- nrow(train); p <- ncol(train)
  if (method=="GEDE"){
    Est <- RobEst(train, ...)
    Xhat <- GEDE(test, Est=Est, predictors=predictors, ...)
  } else if (method=="lasso"){
    Xhat <- matrix(0, nrow(test), p)
    for (i in 1:p) {
      ## use a five-fold CV to select the best lasso model
      Xi.idx <- setdiff(predictors, i)
      mod.i <- cv.glmnet(x=train[,Xi.idx], y=train[,i], type.measure="mse",
                         alpha=1, nfolds=5, ...)
      Xhat[,i] <- predict(mod.i, s=mod.i$lambda.min, newx=test[,Xi.idx])
    }
  } else {
    stop(paste0("The method you specified, ", method, ", has not been implemented in Enhancer() yet."))
  }
  return(Xhat)
}
