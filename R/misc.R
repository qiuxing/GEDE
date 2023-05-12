######################################################################
## misc. useful functions
######################################################################

## total or relative RMSE
rmse <- function(prediction, truth, relative=FALSE) {
  mse <- mean(as.numeric(prediction-truth)^2, na.rm=TRUE)
  if (relative) {
    truth.norm2 <- mean(as.numeric(truth)^2, na.rm=TRUE)
    return(sqrt(mse/truth.norm2))
  } else{ #just return the RMSE
    return(sqrt(mse))
  }
}


## columnwise RMSE, with an option for scaling.
colRMSE <- function(prediction, truth, scaled = FALSE) {
  colrmse <- sqrt(colMeans( (prediction -truth)^2, na.rm=TRUE))
  names(colrmse) <- colnames(truth)
  if (scaled) colrmse <- colrmse/colVars(truth, std=TRUE, na.rm=TRUE)
  return(colrmse)
}

## A convenient function to visualize a matrix
mplot <- function(mat, ncolor=100, col.min=NULL, col.max=NULL, ...) {
  mat <- as.matrix(mat)
  mmin <- min(mat, na.rm=TRUE); mmax <- max(mat, na.rm=TRUE)
  if (is.null(col.min)) col.min <- mmin
  if (is.null(col.max)) col.max <- mmax
  ##
  if (col.min<=mmin & col.max<mmax) {
    breaks <- c(seq(col.min, col.max, length.out=ncolor), mmax)
  } else if (col.min>mmin & col.max<mmax) {
    breaks <- c(mmin, seq(col.min, col.max, length.out=ncolor-1), mmax)
  } else if (col.min>mmin & col.max>=mmax) {
    breaks <- c(mmin, seq(col.min, col.max, length.out=ncolor))
  } else { #col.min<=mmin & col.max>=mmax
    breaks <- seq(col.min, col.max, length.out=ncolor+1)
  }
  image(t(mat[nrow(mat):1,]), col=grey(seq(1, 0, length.out=ncolor)),
        breaks=breaks, xaxt="n", yaxt="n",  bty="n",
        asp=nrow(mat)/ncol(mat),...)
}

## this is a safe replacement for diag(). R's diag() behaves
## differently when x is a vector or a single number
.diag2 <- function(x) diag(x,nrow = length(x))

## A Hampel filter based on Med +/- 3*MAD criterion for Y_{nxp}.
Hampel <- function(Y, nMAD=3, arr.ind=FALSE) {
  if (is.null(nMAD)) { #no filter
    out.idx <- integer(0)
  } else {
    Ymed <- colMedians(Y, na.rm=TRUE)
    Ymad <- colMads(Y, na.rm=TRUE)
    absYc <- abs(sweep(Y, 2, Ymed))
    out.idx <- which(sweep(absYc, 2, nMAD*Ymad)>0, arr.ind=arr.ind)
  }
  return(out.idx)
}

## A robust and fast way to compute beta=(X'X)^{-1} X'Y for highdim Y. We
## assume that the intercept must be included, and we use FWL theorem
## to conduct the regression in two steps. NAs are permitted in Y and
## X. If X==NULL (default, just return colmeans of Y in matrix format.
RobReg <- function(Y, X=NULL, d.prop=1e-6, dmin=1e-9, return.Yhat=FALSE){
  Y <- as.matrix(Y); n <- nrow(Y); m <- ncol(Y)
  Ybars <- colMeans(Y, na.rm=TRUE)
  if (is.null(X)){ #just return the intercepts
    beta0 <- Ybars; betahat <- NULL
  } else { #X is not null
    X <- as.matrix(X); p <- ncol(X)+1
    Xbars <- colMeans(X, na.rm=TRUE)
    ## center Y, and replace NAs in Yc by 0
    Yc <- sweep(Y, 2, Ybars); Yc <- replace(Yc, is.na(Yc), 0)
    ## center X, and replace NAs in Xc by 0
    Xc <- sweep(X, 2, Xbars); Xc <- replace(Xc, is.na(Xc), 0)
    ## Use SVD to conduct numerically robust regression
    o <- svd(Xc)
    ## singular value threshold. Only those eigenvalues that passed this
    ## threshold are used in the hat matrix.
    thresh <- max(sum(o$d) * d.prop, dmin)
    idx <- which(o$d>thresh)
    if (length(idx) < ncol(X)) warning("The design matrix is highly collinear. Some very small singular values used in the regression are omitted.")
    ##
    if (length(idx)==0) { #no singular value is left
      betahat <- matrix(0, p-1, m)
    } else {
      U2 <- o$u[, idx, drop=FALSE]; V2 <- o$v[, idx, drop=FALSE]
      V2D2inv <- sweep(V2, 2, o$d[idx], "/")
      betahat <- V2D2inv%*%(t(U2)%*%Yc)
    }
    rownames(betahat) <- colnames(X); colnames(betahat) <- colnames(Y)
    ## beta0
    beta0 <- Ybars -drop(t(betahat)%*%Xbars)
  }
  ## combine beta0 with betahat
  betahat <- rbind(beta0=beta0, betahat)
  ## compute and return Yhat if return.Yhat=TRUE
  if (return.Yhat) {
    if (is.null(X)) {
      Yhat <- rep(1,n)%*%t(beta0)
    } else {
      Yhat <- cbind(1, X)%*%betahat
    }
    dimnames(Yhat) <- dimnames(Y)
    return(list(betahat=betahat, Yhat=Yhat))
  } else {
    return(betahat)
  }
}

## the predictor for RobReg. newX cannot be NULL.
predict.RobReg <- function(betahat, newX) {
  newX <- as.matrix(newX); n <- nrow(newX)
  if (nrow(betahat)==1) { #only the intercept
    Yhat <- rep(1,n)%*%t(betahat[1,])
  } else {
    Yhat <- cbind(1, newX)%*%betahat
  }
  return(Yhat)
}


## createFolds function copied from caret
createFolds <- function (y, k = 10, list = TRUE, returnTrain = FALSE)
{
  if (inherits(y, "Surv"))
    y <- y[, "time"]
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2)
      cuts <- 2
    if (cuts > 5)
      cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0)
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k,
                                                               size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
                        sep = "")
    if (returnTrain)
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}

## To compute the vector of adjustments for t- and F-tests
t.adj.coef <- function(Est, empirical=TRUE) {
  Tk <- Est$Tk; K <- Est$K; Lk <- Est$Lk
  sigma2 <- Est$sigma2; lks <- Est$lks[1:K]
  Tw1 <- sweep(Tk, 2, sqrt(Lk), "*")
  ## Sigma.jj are variances for the original data
  Sigma.jj <- rowSums(Tw1^2)+sigma2
  ## SigmaTilde.jj are variances for the enhanced data
  if (empirical) {
    Tw2 <- sweep(Tk, 2, Lk*sqrt(lks)/(sigma2+Lk), "*")
  } else {
    Tw2 <- sweep(Tk, 2, Lk/sqrt(sigma2+Lk), "*")
  }
  SigmaTilde.jj <- rowSums(Tw2^2)
  return(sqrt(SigmaTilde.jj/Sigma.jj))
}

## A matrix-valued F-test for comparing two regression models
RegFtest <- function(Y, X, X0=NULL, include.intercept=TRUE, r.adj=1) {
  X <- as.matrix(X); n <- nrow(X)
  if (include.intercept) {
    X <- cbind(1, X)
    X0 <- cbind(rep(1,n), X0)
  }
  ## fit the alternative model
  mod1 <- lm.fit(X, Y); p1 <- mod1$rank
  Yhat1 <- mod1$fitted.values; Res1 <- mod1$residuals
  ## fit the null model
  if (is.null(X0)) { #implies no intercept
    p0 <- 0; Yhat0 <- rep(0, n); Res0 <- Y
  } else {
    X0 <- as.matrix(X0)
    mod0 <- lm.fit(X0, Y); p0 <- mod0$rank
    Yhat0 <- mod0$fitted.values; Res0 <- mod0$residuals
  }
  df1 <- p1-p0; df2 <- mod1$df.residual
  RSS1 <- colSums(Res1^2); RSS0 <- colSums(Res0^2)
  Top <- (RSS0-RSS1)/df1; Bottom <- RSS1/df2
  Bottom.adj <- Bottom/(r.adj^2)
  F=Top/Bottom.adj; pvals <- 1-pf(F, df1, df2)
  rr <- data.frame(RSS0=RSS0, RSS1=RSS1, Top=Top, Bottom=Bottom,
                   Bottom.adj=Bottom.adj, F=F, pvals=pvals)
  rownames(rr) <- colnames(Y)
  ## squeeze df1 and df2 into rr as attributes
  attr(rr, "df1") <- df1; attr(rr, "df2") <- df2
  return(rr)
}

## limma is a wrapper for eBayes. X is a vector or matrix of
## covariates (without the intercept by default). r.adj is the
## adjustment coefficient.
limma <- function(gdata, X, include.intercept=TRUE, r.adj=1){
  Xn0 <- deparse(substitute(X)) #the object name of X
  X <- as.matrix(X)
  if (is.null(colnames(X))) {
    if (ncol(X)==1){
      ## use the object name of the input as variable name
      colnames(X) <- Xn0
    } else {
      colnames(X) <- paste0("X", 1:ncol(X))
    }
  }
  ## remove missing samples
  na.id <- apply(X, 1, function(x) any(is.na(x)))
  gdata <- gdata[, !na.id]; X <- X[!na.id,,drop=FALSE]
  ## combine with the intercepts
  if (include.intercept) X <- cbind(Intercept=1, X)
  ## 
  efit <- eBayes(lmFit(gdata, X))
  betahat <- efit$coefficients
  colnames(betahat) <- paste0("betahat.", colnames(betahat))
  ## the adjusted tstats
  if (include.intercept) {
    tstats=sweep(efit$t[,-1, drop=FALSE], 1, r.adj, "*")
  } else {
    tstats=sweep(efit$t, 1, r.adj, "*")
  }
  Xn <- colnames(tstats)
  colnames(tstats) <- paste0("tstat.", Xn)
  ## p-values computed from tstats
  pvals <- 2*pt(-abs(tstats), df = efit$df.total)
  colnames(pvals) <- paste0("pvals.", Xn)
  adjP <- matrix(p.adjust(pvals, "BH"), nrow=nrow(pvals))
  colnames(adjP) <- paste0("adjP.", Xn)
  return(data.frame(betahat, tstats, pvals, adjP))
}
