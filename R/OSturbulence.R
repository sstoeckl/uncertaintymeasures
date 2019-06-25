#' @title Calculate Out-of-sample turbulence
#'
#' @description Calculate OOS-turbulence based on a cross-section of times-series (xts) data
#'
#'
#' @param X the dataset X in xts format (if not, the program will try to convert it to xts)
#' @param weights optional, can either be a vector (in which case it will be used for each row and must match ncol(X)) or a matrix with the same dimension as X
#' if weights=NULL (the default), equal weights will be generated based on the number of columns (CHECK no of non NA entries per row)
#' @param squared optional, should the square root be taken of the results? (standard: TRUE which correwsponds to squared=FALSE)
#' @param norm optional, should the result be normalized by 1/sum(weights^2) per row (makes expectation equal to 1)
#' @param method optional, if a robust mean/covariance estimator should be used
#' @param s.k optional, lookback window for current observations that are related to the long-term mean
#' @param imp optional, should missing values be imputed? (standard:)
#' @param rolling optional, Should moments be calculated based on a rolling window?
#' @param roll.obs optional, rolling=TRUE: Length of rolling window; rolling=FALSE: Initial estimation window size
#' @param use optional, use="pairwise.complete.obs": What method of dealing with missing values should be used (note, that variables
#' with more than half missing obs are thrown out anyway)
#'
#' @return list containing 4 elements:
#'    - turb: turbulence index
#     - turb.grad: contribution of individual series to turb
#     - mturb: turbulence index based on diagonal covariance matrix
#     - cturb: correlation turbulence turb/mturb
#'
#' @examples
#' require(xts)
#' data(X)
#' turb1 <- OSturbulence(X, s.k=1, rolling=TRUE, roll.obs=36)$turb
#' turb2 <- OSturbulence(X,s.k=12, rolling=TRUE, roll.obs=36)$turb
#' turb3 <- OSturbulence(X,s.k=12, rolling=FALSE, roll.obs=36)$turb
#' plot(turb1, lwd=2, main="(Rolling) Out-of-sample turbulence (standard)")
#' lines(turb2,lwd=2,col="red", main="(Rolling) Turbulence with 12 month lookback period",on=NA)
#' lines(turb3,lwd=2,col="blue", main="(Growing window) Turbulence with 12 month lookback period",on=NA)
#'
#' @import xts
#' @import mice
#' @import fAssets
#' @importFrom zoo rollapplyr
#' @importFrom stats cov
#' @importFrom timeSeries colStats
#'
#' @export
OSturbulence <- function(X, weights=NULL, squared=FALSE, norm=FALSE, method=NULL, s.k=1, imp=FALSE,
                         rolling=FALSE, roll.obs=100, use="complete.obs"){
  # method = c("cov", "mve", "mcd", "MCD", "OGK", "nnve", "shrink", "bagged")
  if (!requireNamespace("xts", quietly = TRUE)) {
    stop("Package \"xts\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("fAssets", quietly = TRUE)) {
    stop("Package \"fAssets\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(imp){
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("Package \"mice\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }
  # check if xts can be produced
  x <- try.xts(X)
  # check for columns that only have NA
  if (any(colSums(!is.na(x))==0)) {stop("Not every column contains at least one non-NA element!")}
  # impute missing values
  if (imp){
    x.cn <- colnames(x); colnames(x) <- NULL
    x <- reclass(complete(mice(x,m = 5, method="pmm", maxit=50, seed=500,printFlag = FALSE)),x)
    colnames(x) <- x.cn
  }
  # initialise output vectors
  x <- if (is.vector(X)) matrix(X, ncol=length(X)) else (as.matrix(X))
  m <- dim(x)[1]; n<-dim(x)[2]
  x.w <- if (is.null(weights)) {sweep(x*0+1,MARGIN = 1,FUN = "/",STATS = rowSums(x*0+1,na.rm=TRUE))} else {as.matrix(weights)}
  x.turb.grad <- x*NA; x.mturb.grad <- x.turb.grad
  x.turb <- x.turb.grad[,1]
  x.mturb <- x.turb
  # dimension warning
  if((roll.obs-1)<=n){stop("You need more observations to start estimation with (cross-section>number of obs)")}
  #
  for (i in roll.obs:m){ #for i>roll.obs we get finally out.of.sample
    # (i-1) for one future observation
    if (rolling==FALSE) {x.sel <- x[1:(i-1),]} else {x.sel <- x[max(1,i-roll.obs):(i-1),]}
    # check if there are cols with more than half NA
    cols <- which(apply(x.sel,2,function(y) sum(!is.na(y)))>=min(n,roll.obs)) # if rolling window is smaller than no of PFs
    if (length(cols)>0){
      if (!is.null(method)) {
        x.moms <- assetsMeanCov(x.sel[,cols],method=method)
        x.Cov <- getCovRob(x.moms)
        x.mean <- getCenterRob(x.moms)
      }	else {
        x.Cov <- (stats::cov(x.sel[,cols],use=use)) ## nearPD removed, replaced by pseudo-inverse
        x.mean <- colMeans(x.sel[,cols],na.rm=TRUE)
      }
      x.iCov <- .mpinv(x.Cov)
      x.dCov <- matrix(0,nrow = nrow(x.Cov),ncol = ncol(x.Cov))
      diag(x.dCov) <- diag(x.Cov)
      x.diCov <- .mpinv(x.dCov)
      # calculate (x*-mu) where x* is the average of the s.k most recent observations
      x.msw <- colStats(x[max(1,(i-s.k+1)):i,cols,drop=FALSE],FUN=mean,align="right",na.pad=TRUE)-x.mean
      # if there are missing row-values, set them to 0 to be able to matrix multiply below
      x.msw[is.na(x[i,cols])] <- 0
      x.w[i,is.na(x[i,cols])] <- 0
      x.w <- sweep(x.w,1,STATS = rowSums(x.w),FUN = "/")#rebalance weights
      #
      x.turb.grad[i,cols] <- ((x.w[i,cols] * x.msw) %*% x.iCov) * (x.w[i,cols] * x.msw)
      x.mturb.grad[i,cols] <- ((x.w[i,cols] * x.msw) %*% x.diCov) * (x.w[i,cols] * x.msw)
      # Sum to get indices
      x.turb[i] <- rowSums(x.turb.grad[i,,drop=FALSE],na.rm=TRUE)
      x.mturb[i] <- rowSums(x.mturb.grad[i,,drop=FALSE],na.rm=TRUE)
      # recover NAs
      x.turb.grad[i,is.na(x[i,cols])] <- NA
      x.mturb.grad[i,is.na(x[i,cols])] <- NA
      x.w[i,is.na(x[i,cols])] <- NA
    } else {next}
  }
  # Normalize
  if (norm==TRUE) {
    x.turb.grad <- sweep(x.turb.grad,MARGIN = 1,FUN = "/",STATS = rowSums(x.w^2,na.rm=TRUE))
    x.turb <- x.turb/rowSums(x.w^2,na.rm=TRUE) # mean of norm(x.turb) is one
    x.mturb <- x.mturb/rowSums(x.w^2,na.rm=TRUE)
  }

  x.cturb <- x.turb/x.mturb
  x.cturb <- reclass(x.cturb,try.xts(X))

  if (squared==FALSE) {
    x.mturb <- x.mturb^0.5 # moved back because we need mturb^2
    x.turb <- x.turb^0.5
    x.turb.grad <- sign(x.turb.grad)*abs(x.turb.grad)^0.5
  }
  x.mturb <- reclass(x.mturb,try.xts(X))
  x.turb <- reclass(x.turb,try.xts(X))
  x.turb.grad <- reclass(x.turb.grad,try.xts(X))
  colnames(x.turb) <- "turb"; colnames(x.mturb) <- "mturb"; colnames(x.cturb) <- "cturb"
  colnames(x.turb.grad) <- colnames(X)
  result <- list(turb =x.turb, turb.grad=x.turb.grad, mturb =x.mturb, cturb =x.cturb)
  return(result)
}
