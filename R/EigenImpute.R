################################################################################
## Below starts the main functions
################################################################################

## 05/01/2023. RobMod() needs to be re-written to deal with covariates


## this function produces regression-like linear predictors of Y given
## X, based on the initial mean/covariance estimator provided by
## RobEst() or a compatible method.  EstObj must be a list that
## contains: K, muhat, Tk, Lk, and sigma2.  Xlist is a *list* of the
## indices of covariates.  It returns a list of linear predictive
## models that can be used in different applications such as
## EigenImpute() and GenePred().
RobMod <- function(EstObj, Xlist) {
  K <- EstObj$K; muhat <- EstObj$muhat
  Tk <- EstObj$Tk; Lk <- EstObj$Lk; sigma2 <- EstObj$sigma2
  mods <- list()
  if (K==0) {
    warning("K=0 means all variables are independent therefore the best estimations are marginal means.")
    for (j in 1:length(Xlist)){
      Xj <- Xlist[[j]]
      beta1.j <- matrix(0, nrow=p-length(Xj), ncol=length(Xj))
      beta0.j <- muhat[-Xj]
      mods[[j]] <- list(beta0=beta0.j, beta1=beta1.j)
    }
  } else if (K==1) {
    ## this is a special case due to two reasons: (a) Tk is a vector;
    ## Lk is a number, so we must use special care to deal with it;
    ## (b) there exists a simple numerical shortcut (see my report)
    ## for this case.
    for (j in 1:length(Xlist)){
      Xj <- Xlist[[j]]
      muhat2 <- muhat[Xj]; muhat1 <- muhat[-Xj]
      Tk2 <- Tk[Xj,,drop=FALSE]; Tk1 <- Tk[-Xj,,drop=FALSE]
      ## d12 is d1^2 in Equation(9)
      d12 <- sum(Tk2^2)
      beta1.j <- Lk/(d12 + sigma2) * (Tk1 %*% t(Tk2))
      beta0.j <- muhat1 - beta1.j %*% muhat2
      mods[[j]] <- list(beta0=beta0.j, beta1=beta1.j)
    }
  } else {
    for (j in 1:length(Xlist)) {
      Xj <- Xlist[[j]]
      muhat2 <- muhat[Xj]; muhat1 <- muhat[-Xj]
      Tk2 <- Tk[Xj,,drop=FALSE]; Tk1 <- Tk[-Xj,,drop=FALSE]
      Tk1b <- sweep(Tk1, 2, sqrt(Lk), "*")
      Tk2b <- sweep(Tk2, 2, sqrt(Lk), "*")
      ss <- svd(Tk2b)
      U <- ss$u; V <- ss$v; dd <- ss$d
      ## beta1.j <- Tk1 %*% .diag2(sqrt(Lk)) %*% V %*% .diag2(dd/(dd^2+sigma2)) %*% t(U)
      beta1.j <- (Tk1b%*%sweep(V, 2, dd/(dd^2+sigma2), "*")) %*% t(U)
      beta0.j <- muhat1 - beta1.j %*% muhat2
      mods[[j]] <- list(beta0=beta0.j, beta1=beta1.j)
    }
  }
  return(mods)
}

############################################################
## several specific applications
############################################################

## EigenImpute() is a multivariate imputation method. Again, note that
## Ymiss is nxp-dimensional, not pxn dimensional as most gene
## expression data are.  Note that this function will not work if some
## samples are completely missing.  It is a good idea to remove
## samples in Ymiss with *too many* missing values before applying this
## function.

## 04/17/2023. A much more memory-efficient algorithm. HD: an
## iterative algorithm that runs faster for high-dim data
EigenImpute <- function(EstObj, Ymiss, predictors=seq(1:ncol(Y)), HD=FALSE, HD.iter=5) {
  ## if no missing, we just return Ymiss
  NA.idx <- which(is.na(Ymiss))
  if (length(NA.idx)==0) {
    warning("There is no missing value in Ymiss. No imputation is required.")
    return(Ymiss)
  }
  ##
  K <- EstObj$K; mumat <- EstObj$mumat
  Tk <- EstObj$Tk; Lk <- EstObj$Lk; sigma2 <- EstObj$sigma2
  n <- nrow(Ymiss); p <- ncol(Ymiss)
  ## gather the information about locations of non-missing data
  Xlist0 <- apply(Ymiss, 1, function(x) which(!is.na(x)))
  ## nx is the number of non-missing variables for each sample
  nx <- sapply(Xlist0, length)
  ## we only need to impute for incomplete samples
  incomplete.cases <- which(nx<p)
  Xlist <- Xlist0[incomplete.cases]
  ## Yc is the centered data to be imputed
  Yc <- Ymiss-mumat
  if (K==0) {
    warning("K=0 means all variables are independent therefore the best estimations are marginal means.")
    Yc[NA.idx] <- 0
  } else { #K >=1; imputation is nontrivial
    if (HD) { #an iterative procedure for high-dim data
      Yc[NA.idx] <- 0
      Tk2 <- Tk[predictors, , drop=FALSE]
      ss <- svd(sweep(Tk2, 2, sqrt(Lk), "*"))
      U <- ss$u; V <- ss$v; dd <- ss$d
      Tktilde <- Tk%*%sweep(V, 1, sqrt(Lk), "*")
      Utilde <- sweep(U, 2, dd/(sigma2+dd^2), "*")
      for (k in 1:HD.iter){
        ## We are over computing (the entire Y2c) here. It could be
        ## done much faster in C.
        Y2c <- Yc[,predictors,drop=FALSE]
        Yc[NA.idx] <- ((Y2c%*%Utilde)%*%t(Tktilde))[NA.idx]
      }
    } else { #the slower but more accurate approach
      for (i in 1:length(incomplete.cases)){
        i0 <- incomplete.cases[i]
        Xi <- intersect(Xlist[[i]], predictors)
        NA.idx.i <- which(is.na(Yc[i0,]))
        y2c <- Yc[i0,Xi] #this is a vector
        Tk2 <- Tk[Xi,,drop=FALSE]; Tk1 <- Tk[NA.idx.i,,drop=FALSE]
        if (K==1) { #a special computational shortcut
          d12 <- sum(Tk2^2)
          Yc[i0, NA.idx.i] <- Lk/(d12 + sigma2)*sum(Tk2*y2c)*Tk1
        } else { #K>=2
          ss <- svd(sweep(Tk2, 2, sqrt(Lk), "*"))
          U <- ss$u; V <- ss$v; dd <- ss$d
          Tk1tilde <- Tk1%*%sweep(V, 1, sqrt(Lk), "*")
          Utilde <- sweep(U, 2, dd/(sigma2+dd^2), "*")
          Yc[i0, NA.idx.i] <- drop(Tk1tilde%*%(t(Utilde)%*%y2c))
        }
      }
    }
  }
  return(Yc+mumat)
}
