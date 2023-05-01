################################################################################
## useful auxiliary functions for the EigenImpute method
################################################################################

## outDetectMD <- function(EstObj, Y) {
##   ## outlier detection, assuming K>1
##   K <- EstObj$K; muhat <- EstObj$muhat
##   Tk <- EstObj$Tk; Lk <- EstObj$Lk; sigma2 <- EstObj$sigma2
##   ResVar <- sapply(1:p, function(i) {
##     Tk.i <- Tk[i,]; Tk.neg.i <- Tk[-i,]
##     beta1.i <- mods[[i]]$beta1
##     drop(Tk.i %*% diag(Lk) %*% ( Tk.i - t(Tk.neg.i)%*% t(beta1.i))) + sigma2
##   })
##   out.z <- sweep(newY - Ypred, 2, sqrt(ResVar), "/")


## }

######################################################################
## Functions related to estimating eigenvalues of the covariance
## matrix
######################################################################

## this function estimates the number of significant PCs.
K.est <- function(lks, l.remain, n, p, method=c("REk", "vprop"), vprop=.8) {
  method <- match.arg(method)
  ## Check if lks are all >=0. Gives warning when not all lks is >=0.
  if (any(lks < 0)) {
    warning("Some eigenvalues are negative. They are replaced by zeros.")
    lks <- replace(lks, lks<0, 0)
  }
  ## sort lks
  lks <- sort(lks, decreasing=TRUE)
  tvar <- sum(lks)+l.remain #total variance
  if (method=="vprop"){
    vp <- cumsum(lks)/tvar
    mvp <- vp[length(vp)] #max varprop
    if (mvp<vprop) {
      warning(paste0("The maximum proportion of variance explained is max(varprop)=", round(mvp, 3), ", which is less than the given vprop cutoff, vprop=", vprop, ". As a result, Kstar is set to be the length of lks."))
      Kstar <- length(lks)
    } else {
      Kstar <- which(vp>=vprop)[1]
    }
  } else if (method=="REk") {
    pstar <- min(p, n-1)
    ## Now compute the REk stats
    tk <- sapply(0:pstar, function(k) p*((p-k)*sum( (lks[-(1:k)])^2) / sum(lks[-(1:k)])^2 -(1+p/n)) - p/n)
    REk <- sapply(0:pstar, function(k) 1/4*(n/p)^2 * tk[k+1]^2 + 2*(k+1))
    Kstar <- which.min(REk)-1
  }
  return(Kstar)
}

## We assume that Y contains no missing value nor outliers. px: DF of
## the covariates
SimpleEst <- function(Y, y.centered=TRUE, px=1, K="auto", K.method=c("REk", "vprop"), vprop=0.8, Kmax=100) {
  n <- nrow(Y); p <- ncol(Y); K.method <- match.arg(K.method)
  pstar <- min(p, n-px)
  if (y.centered) { #Y is already centered
    mumat <- matrix(0, n, p); Yc <- Y
  } else { #we need to center Y
    mumat <- rep(1,n)%*%t(colMeans(Y))
    Yc <- Y-mumat
  }
  ## ss <- svd(Yc)
  ## Tmat <- ss$v; lks <- ss$d^2/(n-1)
  ss <- hd.eigen(Yc, center=FALSE, scale=FALSE, k=min(Kmax,pstar), vectors=TRUE, large=TRUE)
  Tmat <- ss$vectors; lks <- ss$values
  ## l.remain is the variance unexplained by the Kmax eigenvalues
  l.remain <- round(sum(Yc^2)/(n-px)-sum(lks), 4)
  ## for convenience, output varprops for screeplot
  varprops <- cumsum(lks)/(sum(lks)+l.remain)
  if (K=="auto"){
    K=K.est(lks, l.remain, n=n, p=p, method=K.method, vprop=vprop)
  } else if( K<0) {
    stop("K must be a non-negative integer or auto.")
  }
  ## Now we have a nonnegative K
  if (K==0) {
    warning("Either you or the automatic model selection procedure selected K=0 or , which implies that all variables are independent. No information can be borrowed from other variables, therefore the best prediction is mean imputation. If that's not expected, consider manually select a reasonable K instead.")
    Tk <- Lk <- NA; sigma2 <- mean(lks); PCs <- NA
  } else { #K >= 1
    if (K>=pstar) {
      warning("Your K is greater or equal to min(p,n-1), which implies that sigma2=0.  Consider using a smaller K.")
      K <- pstar; Tk <- Tmat; Lk <- lks; sigma2 <- 0
    } else { # 0<K<pstar; the main case
      Tk <- Tmat[, 1:K, drop=FALSE]
      ## estimate the first K eigenvalues
      sigma2 <- pstar/(p*(pstar-K)) *sum(lks[-(1:K)])
      Lk <- lks[1:K]-p/pstar*sigma2
    }
    PCs <- Yc%*%Tk
    rownames(PCs) <- rownames(Y); colnames(PCs) <- paste0("PC", 1:K)
  }
  return(list(mumat=mumat, lks=lks, l.remain=l.remain, varprops=varprops, Tk=Tk, K=K, Lk=Lk, sigma2=sigma2, PCs=PCs))
}

## covariates must be numeric
RobEst <- function(Y, covariates=NULL, K="auto", K.method=c("REk", "vprop"), vprop=0.8, Kmax=100, nMAD=3, HD=FALSE, HD.iter=5) {
  n <- nrow(Y); p <- ncol(Y)
  ## 1. Outlier removal
  out.idx <- Hampel(Y, nMAD=nMAD)
  Ymiss <- Y; Ymiss[out.idx] <- NA
  ## 2. Estimate the mean values by a robust and fast regression function
  mumat <- RobReg(Y, covariates)
  ## 3. Initial parameter estimation based on Yc0
  Yc0 <- Ymiss-mumat
  Yc0 <- replace(Yc0, is.na(Yc0), 0) #replace NAs by 0
  K.method <- match.arg(K.method)
  suppressWarnings( Est0 <- SimpleEst(Yc0, K=K, K.method=K.method, vprop=vprop, Kmax=Kmax) )
  Est0$mumat <- mumat #manually add mumat back to Est0
  ## 4. Impute missing values 
  Y1 <- EigenImpute(Est0, Ymiss, HD=HD, HD.iter=HD.iter)
  mumat <- RobReg(Y1, covariates)
  Y1c <- Y1-mumat
  ## 5. Second (final) round of parameter estimation
  Est1 <- SimpleEst(Y1c, K=K, K.method=K.method, vprop=vprop, Kmax=Kmax)
  ## 5. append some useful information to Est1 before return
  Est1$mumat <- mumat
  Est1[["NA.idx"]] <- which(is.na(Y)); Est1[["out.idx"]] <- out.idx
  return(Est1)
}
