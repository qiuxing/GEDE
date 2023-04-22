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
K.est <- function(lks, n, p="auto", tol=1e-4){
  ## if p="auto", assume that length of lks is p
  if (p=="auto"){
    p <- length(lks)
  } else {
    ## otherwise, append zeros to the rest of lks
    if (p> length(lks)) lks <- c(lks, rep(0, p-length(lks)))
  }
  ## Check if lks is >=0. Gives warning when not all lks is >=0.
  if (any(lks < 0)){
    warning("Some eigenvalues are negative. They are replaced by zeros.")
    lks <- ifelse(lks<0, 0, lks)
  }
  ## Also warns if number of >0 eigenvalues is < min(n,p)
  N <- min(n,p)-1
  if (sum(lks>0) < N) {
    warning(paste0("Number of positive eigenvalues is ", sum(lks>0), ", which is less than min(n,p)-1 = ", N, ". I have to use ad hoc workaround to obtain AIC and MDL values, which may not be accurate! "))
  }
  ## sort lks
  lks <- sort(lks, decreasing=TRUE)
  ## An ad hoc upper bound. We will not consider eigenvalues with very
  ## small variance.
  tol <- sum(lks)*tol
  N2 <- ifelse(min(lks)<tol, which(lks<tol)[1], N)
  ## Now computes some important summary stats
  gk <- sapply(0:N2, function(k) exp(mean(log(lks[-(1:k)]))))
  ak <- sapply(0:N2, function(k) mean(lks[-(1:k)]))
  tk <- sapply(0:N2, function(k) p*((p-k)*sum( (lks[-(1:k)])^2) / sum(lks[-(1:k)])^2 -(1+p/n)) - p/n)
  AICk <- sapply(0:N2, function(k) -2*(p-k)*n*log(gk[k+1]/ak[k+1]) + 2*k*(2*p-k))
  MDLk <- sapply(0:N2, function(k) -(p-k)*n*log(gk[k+1]/ak[k+1]) + 1/2*k*(2*p-k)*log(n))
  REk <- sapply(0:N2, function(k) 1/4*(n/p)^2 * tk[k+1]^2 + 2*(k+1))
  stats <- cbind(AICk=AICk, MDLk=MDLk, REk=REk)
  rownames(stats) <- paste0("K=", 0:N2)
  return(list(Kstar=c(AICk=which.min(AICk)-1, MDLk=which.min(MDLk)-1, REk=which.min(REk)-1), statistics=stats)
         )
}

## ## We assume that Y contains no missing value nor outliers
SimpleEst <- function(Y, K="auto", Kmax=100) {
  n <- nrow(Y); p <- ncol(Y)
  pstar <- min(p, n-1)
  muhat <- colMeans(Y); Yc <- sweep(Y, 2, muhat)
  ## ss <- svd(Yc)
  ## Tmat <- ss$v; lks <- ss$d^2/(n-1)
  ss <- hd.eigen(Yc, center=FALSE, scale=FALSE, k=min(Kmax,pstar), vectors=TRUE, large=TRUE)
  Tmat <- ss$vectors; lks <- ss$values
  ## l.remain is the variance unexplained by the Kmax eigenvalues
  l.remain <- round(sum(Yc^2)/(n-1)-sum(lks), 4)
  ## for convenience, output varprops for screeplot
  varprops <- cumsum(lks)/(sum(lks)+l.remain)
  if (K=="auto"){
    K=K.est(lks, n=n, p=p)[["Kstar"]][["REk"]]
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
  return(list(muhat=muhat, lks=lks, l.remain=l.remain, varprops=varprops, Tk=Tk, K=K, Lk=Lk, sigma2=sigma2, PCs=PCs))
}

## 04/14/2023. 
RobEst <- function(Y, K="auto", Kmax=100, nMAD=3, HD=FALSE, HD.iter=5) {
  ## 1. Outlier removal
  out.idx <- Hampel(Y, nMAD=nMAD)
  Ymiss <- Y; Ymiss[out.idx] <- NA
  ## 2. Initial parameter estimation based on Yc0
  Ybar <- colMeans(Ymiss, na.rm=TRUE)
  Yc0 <- sweep(Ymiss, 2, Ybar)
  Yc0 <- replace(Yc0, is.na(Yc0), 0) #replace NA by 0
  Est0 <- SimpleEst(Yc0, K=K, Kmax=Kmax)
  Est0$muhat <- Ybar #manually add Ybar back to Est0
  ## 3. 
  Y1 <- EigenImpute(Est0, Ymiss, HD=HD, HD.iter=HD.iter)
  ## 4. Second (final) round of parameter estimation
  Est1 <- SimpleEst(Y1, K=K, Kmax=Kmax)
  ## 5. append some useful information to Est1 before return
  Est1[["NA.idx"]] <- which(is.na(Y)); Est1[["out.idx"]] <- out.idx
  return(Est1)
}