#' @title Calculate (in-sample) the coefficient to adjust the optimal Markowitz portfolio rule based on financial turbulence
#'
#' @description Adjust the optimal Markowitz Portfolio rule according to model 1.2.2 of Garlappi, Uppal & Wang (2008, RFS)
#' According tom the model, we use the minimax rule of Gilboa & Schmeidler (1989) to use the worst (min/inf) of all possible
#' priors for the calculation of the optimal portfolio rule.
#' Uncertainty about expected returns estimated jointly for all assets
#' This coefficient explains how much to move out of the provided factor towards the risk-free rate or (if using the model without the risk-free rate)
#' towards the global Minimum-Variance portfolio
#' The measure of current uncertainty is based on Stöckl (2019, WP)
#'
#'
#' @param X the dataset X in xts format (if not, the program will try to convert it to xts). Monthly returns.
#' @param weights optional, can either be a vector (in which case it will be used for each row and must match ncol(X)) or a matrix with the same dimension as X
#' if weights=NULL (the default), equal weights will be generated based on the number of columns (CHECK no of non NA entries per row)
#' @param method optional, if a robust mean/covariance estimator should be used
#' @param s.k optional, lookback window for current observations that are related to the long-term mean
#' @param imp optional, should missing values be imputed? (standard: FALSE)
#' @param rf optional, should the case with (default) or without a risk-free rate be used (the latter case is not implemented)
#'
#' @return timeseries/vector
#     - sturb: how much should be invested in the provided index, given its Sharpe-ratio, turbulence (paremeter uncertainty) and investors risk-aversion
#'
#' @examples
#' require(xts)
#' data(X)
#' sturb1 <- ISminimax_one(X, s.k=1)
#' sturb2 <- ISminimax_one(X,s.k=12)
#' sturb3 <- ISminimax_one(X,s.k=12, method="MCD")
#' plot(sturb1$sturb, lwd=2, main="Portfolio Adjustment coefficient with risk-free rate",ylim=c(-0.1,1))
#' lines(sturb2$sturb,lwd=2,col="red", main="Portfolio Adjustment coefficient based on turbulence with 12 month lookback period",on=NA,ylim=c(-0.1,1))
#' lines(sturb3$sturb,lwd=2,col="blue", main="Portfolio Adjustment coefficient based on turbulence with 12 month lookback period and robust MCD estimator",on=NA,ylim=c(-0.1,1))
#'
#' @import xts
#' @import mice
#' @import fAssets
#' @importFrom zoo rollapplyr
#' @importFrom stats cov
#' @importFrom timeSeries colStats
#'
#' @export
ISminimax_one <- function(X, weights=NULL, method=NULL, s.k=1, imp=FALSE, rf=TRUE){
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
  # Calculate full-sample turbulence
  x <- if (is.vector(X)) matrix(X, ncol=length(X)) else (as.matrix(X))
  m <- dim(x)[1]; n<-dim(x)[2]
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
  # now create the (x-mu)
  #x.msw <- sweep(x,MARGIN = 2,FUN = "-",STATS = x.mean)
  # calculate (x*-mu) where x* is the average of the s.k most recent observations
  x.msw <- sweep(rollapplyr(x,by.column = TRUE,width = s.k, FUN=mean, fill = NA),MARGIN = 2,FUN = "-",STATS = x.mean)
  # if there are missing row-values, set them to 0 to be able to matrix multiply below
  x.msw[is.na(x.msw)] <- 0
  # Full Sample Sharpe ratio
  x.SR <- x.mean %*% x.iCov %*% x.mean
  # Formula (56) in Kan & Zhou (2007)
  # Subsample version (always use F-correction despite it being wrong for s.k=1)
  x.sturb <- x.msw[,1]*NA
  if (rf){
    if (s.k==1){
      x.epsilon <- ( (m-1)*(m-n+1)/(m-n)/m^2 * rowSums((x.msw %*% x.iCov) * x.msw))
      x.sturb <- apply(as.matrix(1-(x.epsilon/c(x.SR))^0.5),1,FUN=function(x){max(x,0)})
    } else {
      x.epsilon <- ( (m-1)*s.k*(m-m/s.k-n+1)/(s.k-1)/(m/s.k-1)/(m-n)/m * rowSums((x.msw %*% x.iCov) * x.msw))
      x.epsilon[1:(s.k-1)] <- NA
      x.sturb <- apply(as.matrix(1-(x.epsilon/c(x.SR))^0.5),1,FUN=function(x){max(x,0)})
    }
  } else {
    stop("the version without risk-free rate is not yet implemented!\n",)
  }
  x.sturb <- reclass(x.sturb,try.xts(X))
  x.epsilon <- reclass(x.epsilon,try.xts(X))
  return(list(sturb=x.sturb,epsilon=x.epsilon))
}
