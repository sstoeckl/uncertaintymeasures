#' @title Calculate In-sample turbulence
#'
#' @description Calculate IS-turbulence based on a cross-section of times-series (xts) data
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
#' turb1 <- ISturbulence(X, s.k=1)$turb
#' turb2 <- ISturbulence(X,s.k=12)$turb
#' turb3 <- ISturbulence(X,s.k=12, method="MCD")$turb
#' plot(turb1, lwd=2, main="In-sample turbulence (standard)")
#' lines(turb2,lwd=2,col="red", main="Turbulence with 12 month lookback period",on=NA)
#' lines(turb3,lwd=2,col="blue", main="Turbulence with 12 month lookback period and robust MCD estimator",on=NA)
#'
#' @import xts
#' @import mice
#' @import fAssets
#' @importFrom zoo rollapplyr
#' @importFrom stats cov
#' @importFrom timeSeries colStats
#'
#' @export
ISturbulence <- function(X, weights=NULL, squared=FALSE, norm=FALSE, method=NULL, s.k=1, imp=FALSE){
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
  x.turb <- NULL
  x.mturb <- NULL
  #
  if (!is.null(method)) {
    x.moms <- assetsMeanCov(x,method=method)
    x.Cov <- getCovRob(x.moms)
    x.mean <- getCenterRob(x.moms)
  }	else {
    x.Cov <- (stats::cov(x,use="complete")) ## nearPD removed, replaced by pseudo-inverse
    x.mean <- colMeans(x,na.rm=TRUE)
  }
  x.iCov <- .mpinv(x.Cov)
  x.dCov <- 0*x.Cov
  diag(x.dCov) <- diag(x.Cov)
  x.diCov <- .mpinv(x.dCov)
  # now create the (x-mu)
  #x.msw <- sweep(x,MARGIN = 2,FUN = "-",STATS = x.mean)
  # calculate (x*-mu) where x* is the average of the s.k most recent observations
  x.msw <- sweep(rollapplyr(x,by.column = TRUE,width = s.k, FUN=mean, fill = NA),MARGIN = 2,FUN = "-",STATS = x.mean)
  # if there are missing row-values, set them to 0 to be able to matrix multiply below
  x.msw[is.na(x.msw)] <- 0
  x.w[is.na(x.w)] <- 0
  x.w <- sweep(x.w,1,STATS = rowSums(x.w),FUN = "/")#rebalance weights
  x.turb.grad <- ((x.w * x.msw) %*% x.iCov) * (x.w * x.msw)
  x.mturb.grad <- ((x.w * x.msw) %*% x.diCov) * (x.w * x.msw)
  # Sum to get indices
  x.turb <- rowSums(x.turb.grad)
  x.mturb <- rowSums(x.mturb.grad)
  # recover NAs
  x.turb.grad[is.na(x)] <- NA
  x.mturb.grad[is.na(x)] <- NA
  x.w[is.na(x)] <- NA
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
    x.turb.grad <- sweep(x.turb.grad,MARGIN = 1,FUN = "/",STATS = x.turb)  #sign(x.turb.grad)*abs(x.turb.grad)^0.5
  }
  x.mturb <- reclass(x.mturb,try.xts(X))
  x.turb <- reclass(x.turb,try.xts(X))
  x.turb.grad <- reclass(x.turb.grad,try.xts(X))
  colnames(x.turb) <- "turb"; colnames(x.mturb) <- "mturb"; colnames(x.cturb) <- "cturb"
  colnames(x.turb.grad) <- colnames(X)
  result <- list(turb =x.turb, turb.grad=x.turb.grad, mturb =x.mturb, cturb =x.cturb)
  return(result)
}
