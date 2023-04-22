######################################################################
## misc. useful functions
######################################################################

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


