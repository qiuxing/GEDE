################################################################################
## useful auxiliary functions for the EigenImpute method
################################################################################

## this is a safe replacement for diag(). R's diag() behaves
## differently when x is a vector or a single number
.diag2 <- function(x) diag(x,nrow = length(x))

## El1 when sigma2 == 1.
.efun <- function(n,p) {
  ((sqrt(n-2)+sqrt(p))^2 - 1.2065*(sqrt(n-2)+sqrt(p))*(1/sqrt(n-2)+1/sqrt(p))^(1/3))/(n-1)
}


## the expectation of l_1 when sigma^2_epsilon==1
.cfun <- function(n,p,k){
  pstar <- min(p,n-1); enp <- .efun(n,p)
  return(enp - (k-1)/(pstar-1)*(2*enp-2*p/pstar))
}

## the *tail* sum of Elk
.sfun <- function(n,p,K){
  pstar <- min(p,n-1); enp <- .efun(n,p)
  return(p - K*enp + K*(K-1)/(pstar*(pstar-1)) *(enp*pstar-p))
}


######################################################################
## Functions related to estimating eigenvalues of the covariance
## matrix
######################################################################

## this function estimates the number of significant PCs.
K.est <- function(LL, n, p="auto", tol=1e-3){
  ## if p="auto", assume that length of LL is p
  if (p=="auto"){
    p <- length(LL)
  } else {
    ## otherwise, append zeros to the rest of LL
    if (p> length(LL)) LL <- c(LL, rep(0, p-length(LL)))
  }
  ## Check if LL is >=0. Gives warning when not all LL is >=0.
  if (any(LL < 0)){
    warning("Some eigenvalues are negative. They are replaced by zeros.")
    LL <- ifelse(LL<0, 0, LL)
  }
  ## Also warns if number of >0 eigenvalues is < min(n,p)
  N <- min(n,p)-1
  if (sum(LL>0) < N) {
    warning(paste0("Number of positive eigenvalues is ", sum(LL>0), ", which is less than min(n,p)-1 = ", N, ". I have to use ad hoc workaround to obtain AIC and MDL values, which may not be accurate! "))
  }
  ## sort LL
  LL <- sort(LL, decreasing=TRUE)
  ## An ad hoc upper bound. We will not consider eigenvalues with very
  ## small variance.
  tol <- sum(LL)*tol
  N2 <- ifelse(min(LL)<tol, which(LL<tol)[1]-2, N)
  ## Now computes some important summary stats
  gk <- sapply(0:N2, function(k) exp(mean(log(LL[-(1:k)]))))
  ak <- sapply(0:N2, function(k) mean(LL[-(1:k)]))
  tk <- sapply(0:N2, function(k) p*((p-k)*sum( (LL[-(1:k)])^2) / sum(LL[-(1:k)])^2 -(1+p/n)) - p/n)
  AICk <- sapply(0:N2, function(k) -2*(p-k)*n*log(gk[k+1]/ak[k+1]) + 2*k*(2*p-k))
  MDLk <- sapply(0:N2, function(k) -(p-k)*n*log(gk[k+1]/ak[k+1]) + 1/2*k*(2*p-k)*log(n))
  REk <- sapply(0:N2, function(k) 1/4*(n/p)^2 * tk[k+1]^2 + 2*(k+1))
  stats <- cbind(AICk=AICk, MDLk=MDLk, REk=REk)
  rownames(stats) <- paste0("K=", 0:N2)
  return(list(Kstar=c(AICk=which.min(AICk)-1, MDLk=which.min(MDLk)-1, REk=which.min(REk)-1), statistics=stats)
         )
}


## eigenvalue estimators. Right now the "simple" method works better
## than "LBCE".  I plan to improve LBCE in the near future.
LambdaEst <- function(lks, n, p, K, method=c("simple", "LBCE")) {
  method <- match.arg(method)
  pstar <- min(p, n-1)
  if (method=="simple"){
    sigma2hat <- pstar/(p*(pstar-K)) *sum(lks[-(1:K)])
    lambdahat <- lks[1:K]-p/pstar*sigma2hat
  } else if (method=="LBCE") {
    sigma2hat <- sum(lks[-(1:K)])/ .sfun(n,p,K)
    cc <- sapply(1:K, function(k) .cfun(n,p,k))
    lambdahat <- lks[1:K]-sigma2hat*cc
  }
  return(list(lambdahat=lambdahat, sigma2hat=sigma2hat))
}


## InitialEst() is responsible for estimating \mu, the marginal means;
## and \Sigma^(0), a rough initial estimate of covariance matrix.
## mvnmle (the mlest() function) does not work very well. It cannot
## handle high-throughput (p>=50) data and has severe numerical
## issues. Right now it defaults to the basic estimator.
InitialEst <- function(Y, method=c("basic", "mvnmle", "softImpute"), ...) {
  method <- match.arg(method)
  if (method=="basic") {
    muhat <- colMeans(Y, na.rm=TRUE)
    Sigma0 <- cov(Y, use="pairwise.complete.obs")
  } else if (method=="mvnmle") {
    mod <- mlest(Y, ...)
    muhat <- mod$muhat
    Sigma0 <- mod$sigmahat
  } else if (method=="softImpute") {
    stop ("Unavailable at this moment.")
  } else {
    stop ("Only basic, mvnmle, and softImpute are currently implemented.")
  }
  return(list(muhat=muhat, Sigma0=Sigma0))
}

## This function takes Sigma0, and re-estimate / decompose it based on
## the PPCA model.  It returns the eigen-vector matrix (Tk),
## eigenvalues (Lk), and the variance of iid noise (sigma2).  Note
## that Y is supposed to be an nxp-dimensional matrix, so a common
## gene expression matrix has to be transposed to form Y.
RobEst <- function(Y, method1=c("basic", "mvnmle", "softImpute"), method2=c("simple", "LBCE"), K="auto", ...) {
  n <- nrow(Y); p <- ncol(Y); p2 <- min(p, n-1)
  method1 <- match.arg(method1); method2 <- match.arg(method2)
  ## 1. Initial est.
  est1 <- InitialEst(Y, method=method1)
  muhat <- est1$muhat; Sigma0 <- est1$Sigma0
  ## 2. decomposition of Sigma0
  ee <- eigen(Sigma0, symmetric=TRUE); lks <- ee$values; Tmat <- ee$vectors
  if (K=="auto"){
    K=K.est(lks, n=n, p=p)[["Kstar"]][["REk"]]
  } else if( K<0) {
    stop("K must be a non-negative integer or auto.")
  }
  ## Now we have a nonnegative K
  if (K==0) {
    warning("Either you or the automatic model selection procedure selected K=0 or , which implies that all variables are independent. No information can be borrowed from other variables, therefore the best prediction is mean imputation. If that's not expected, consider manually select a reasonable K instead.")
    Tk <- Lk <- NA; sigma2 <- mean(lks)
  } else if (K==p) {
    Tk <- Tmat; Lk <- lks; sigma2 <- 0
  } else { # 0<K<p; the main case
    Tk <- Tmat[, 1:K, drop=FALSE]
    ## estimate the first K eigenvalues
    lk.est <- LambdaEst(lks, n, p, K, method=method2)
    Lk <- lk.est$lambdahat; sigma2 <- lk.est$sigma2hat
  }
  return(list(muhat=muhat, Tmat=Tmat, lks=lks, Tk=Tk, K=K, Lk=Lk, sigma2=sigma2))
}

