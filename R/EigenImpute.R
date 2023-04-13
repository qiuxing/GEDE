################################################################################
## Below starts the main functions
################################################################################

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
      ss <- svd(Tk2 %*% .diag2(sqrt(Lk)))
      U <- ss$u; V <- ss$v; dd <- ss$d
      beta1.j <- Tk1 %*% .diag2(sqrt(Lk)) %*% V %*% .diag2(dd/(dd^2+sigma2)) %*% t(U)
      beta0.j <- muhat1 - beta1.j %*% muhat2
      mods[[j]] <- list(beta0=beta0.j, beta1=beta1.j)
    }
  }
  return(mods)
}

############################################################
## several specific applications
############################################################

## a multivariate imputation method. Again, note that Ymiss is
## nxp-dimensional, not pxn dimensional as most gene expression data
## are.  Note that this function will not work if some samples are
## completely missing.  It is a good idea to remove samples in Ymiss
## with too many missing values before applying this function
EigenImpute <- function(EstObj, Ymiss) {
  n <- nrow(Ymiss); p <- ncol(Ymiss)
  ## gather the information about locations of non-missing data
  Xlist0 <- apply(Ymiss, 1, function(x) which(!is.na(x)))
  ## nx is the number of non-missing variables for each sample
  nx <- sapply(Xlist0, length)
  ## we only need to impute for incomplete samples
  incomplete.cases <- which(nx<p)
  Xlist <- Xlist0[incomplete.cases]
  ## build predictive models
  mods <- RobMod(EstObj, Xlist)
  ## use these linear models to predict missing values
  Y2 <- Ymiss
  for (i in 1:length(incomplete.cases)){
    i0 <- incomplete.cases[i]
    Yi <- Ymiss[i0,]; Xi <- Yi[Xlist[[i]]]
    NA.idx.i <- which(is.na(Yi))
    beta0.i <- mods[[i]]$beta0; beta1.i <- mods[[i]]$beta1
    Yi.hat <- drop(beta0.i + beta1.i %*% Xi)
    Y2[i0, NA.idx.i] <- Yi.hat
  }
  return(Y2)
}

## a way to predict every gene from other genes. Number of covariates
## in the training data Y and newY must be the same (sample sizes can
## be different)
GeneGenePred.old <- function(EstObj, newY) {
  p <- ncol(newY)
  ## for each gene i, we build a predictive model with all other genes
  Xlist <- lapply(1:p, function(i) seq(1,p)[-i])
  mods <- RobMod(EstObj, Xlist)
  ## now apply those linear predictive models to newY
  Ypred <- sapply(1:p, function(i) {
    beta0.i <- mods[[i]]$beta0; beta1.i <- mods[[i]]$beta1
    drop(beta0.i) + drop(newY[,-i,drop=FALSE] %*% t(beta1.i))
  })
  return(Ypred)
}

## 04/01/2020. Added an outlier detection algorithm in this function
GeneGenePred <- function(EstObj, newY, out.detect=FALSE) {
  p <- ncol(newY)
  ## for each gene i, we build a predictive model with all other genes
  Xlist <- lapply(1:p, function(i) seq(1,p)[-i])
  mods <- RobMod(EstObj, Xlist)
  ## now apply those linear predictive models to newY
  Ypred <- sapply(1:p, function(i) {
    beta0.i <- mods[[i]]$beta0; beta1.i <- mods[[i]]$beta1
    drop(beta0.i) + drop(newY[,-i,drop=FALSE] %*% t(beta1.i))
  })
  if (out.detect) {
    ## outlier detection, assuming K>1
    K <- EstObj$K; muhat <- EstObj$muhat
    Tk <- EstObj$Tk; Lk <- EstObj$Lk; sigma2 <- EstObj$sigma2
    if (K < 2){
      stop("I have not implemented K=0 and K=1 yet.")
    } else {
      ResVar <- sapply(1:p, function(i) {
        Tk.i <- Tk[i,]; Tk.neg.i <- Tk[-i,]
        beta1.i <- mods[[i]]$beta1
        drop(Tk.i %*% diag(Lk) %*% ( Tk.i - t(Tk.neg.i)%*% t(beta1.i))) + sigma2
      })
    }
    out.z <- sweep(newY - Ypred, 2, sqrt(ResVar), "/")
      return(list(Ypred=Ypred, out.z=out.z, ResVar=ResVar))
  } else {
    return(Ypred)
  }
}

## A fast implementation specifically designed for gene-gene
## prediction




## more applications!!
