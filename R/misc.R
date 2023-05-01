######################################################################
## misc. useful functions
######################################################################

## total RMSE
rmse <- function(prediction, truth) sqrt(mean(as.numeric(prediction-truth)^2, na.rm=TRUE))


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
  ## Ymed <- apply(Y, 2, median, na.rm=TRUE)
  ## Ymad <- apply(Y, 2, mad, na.rm=TRUE)
  Ymed <- colMedians(Y, na.rm=TRUE)
  Ymad <- colMads(Y, na.rm=TRUE)
  absYc <- abs(sweep(Y, 2, Ymed))
  out.idx <- which(sweep(absYc, 2, nMAD*Ymad)>0, arr.ind=arr.ind)
  return(out.idx)
}

## A robust and fast way to compute X (X'X)^{-1} X'Y for highdim Y. We
## assume that the intercept must be included, and we use FWL theorem
## to conduct the regression in two steps. NAs are permitted in Y and
## X. If X==NULL (default, just return colmeans of Y in matrix format.
RobReg <- function(Y, X=NULL, d.prop=1e-6, dmin=1e-9){
  Y <- as.matrix(Y); n <- nrow(Y)
  Ybars <- colMeans(Y, na.rm=TRUE)
  if (is.null(X)){
    Yhat <- rep(1,n)%*%t(Ybars)
  } else { #X is not null
    X <- as.matrix(X)
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
      Ychat <- matrix(0, n, ncol(Y))
    } else {
      U <- o$u[, idx, drop=FALSE]
      Ychat <- U%*%(t(U)%*%Yc)
    }
    ## Add Ybars back
    Yhat <- sweep(Ychat, 2, Ybars, "+")
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

## limma is a wrapper for eBayes. "v" is a matrix of covariates
## WITHOUT the intercept.
limma <- function(gdata, v){
  v <- as.matrix(v); vn <- colnames(v)
  if (is.null(vn)) {
    vn <- paste0("X", 1:ncol(v)); colnames(v) <- vn
  }
  ## remove missing samples
  na.id <- apply(v, 1, function(x) any(is.na(x)))
  gdata <- gdata[, !na.id]; v <- v[!na.id,,drop=FALSE]
  design.mat <- cbind(Intercept=1, v)
  efit <- eBayes(lmFit(gdata, design.mat))
  betahat <- efit$coefficients
  colnames(betahat) <- paste0("betahat.", colnames(betahat))
  tstats=efit$t[,-1, drop=FALSE]
  colnames(tstats) <- paste0("tstat.", colnames(tstats))
  pvals <- efit$p.value[,-1, drop=FALSE]
  colnames(pvals) <- paste0("pvals.", vn)
  adjP <- matrix(p.adjust(pvals, "BH"), nrow=nrow(pvals))
  colnames(adjP) <- paste0("adjP.", vn)
  return(data.frame(betahat, tstats, pvals, adjP))
}
