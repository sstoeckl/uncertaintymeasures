#' @title Calculate (out-of-sample) the coefficient to adjust the optimal Markowitz portfolio rule based on financial turbulence
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
#' @param X the dataset X in xts format (if not, the program will try to convert it to xts)
#' @param weights optional, can either be a vector (in which case it will be used for each row and must match ncol(X)) or a matrix with the same dimension as X
#' if weights=NULL (the default), equal weights will be generated based on the number of columns (CHECK no of non NA entries per row)
#' @param method optional, if a robust mean/covariance estimator should be used c("cov", "mve", "mcd", "MCD", "OGK", "nnve", "shrink", "bagged")
#' @param s.k optional, lookback window for current observations that are related to the long-term mean
#' @param imp optional, should missing values be imputed? (standard: FALSE)
#' @param rolling optional, Should moments be calculated based on a rolling window?
#' @param roll.obs optional, rolling=TRUE: Length of rolling window; rolling=FALSE: Initial estimation window size
#' @param rf optional, should the case with (default) or without a risk-free rate be used (the latter case is not implemented)
#' @param NOUT optional, if the epsilons are needed as input for a multi-uncertainty portfolio, supply the number of assets per uncertainty group here
#'
#' @return xts/vector
#'     - sturb: how much should be invested in the provided index, given its Sharpe-ratio, turbulence (paremeter uncertainty) and investors risk-aversion
#'     - epsilon: Level of parameter uncertainty
#'
#' @examples
#' require(xts)
#' data(X)
#' sturb1 <- OSminimax_one(X, s.k=1, rolling=TRUE, roll.obs=36)$sturb
#' sturb2 <- OSminimax_one(X,s.k=12, rolling=TRUE, roll.obs=36)$sturb
#' sturb3 <- OSminimax_one(X,s.k=12, rolling=FALSE, roll.obs=36)$sturb
#' plot(sturb1, lwd=2, main="Portfolio Adjustment coefficient with risk-free rate",ylim=c(-0.1,1))
#' lines(sturb2,lwd=2,col="red", main="Portfolio Adjustment coefficient based on turbulence with 12 month lookback period",on=NA,ylim=c(-0.1,1))
#' lines(sturb3,lwd=2,col="blue", main="Portfolio Adjustment coefficient based on turbulence with 12 month lookback period and robust MCD estimator",on=NA,ylim=c(-0.1,1))
#' # relation between in-sample and out-of-sample
#' isturb <- ISminimax_one(X, s.k=12)
#' osturb <- OSminimax_one(X, s.k=12, rolling=FALSE, roll.obs=100)
#' plot(cbind(isturb$sturb,osturb$sturb),col=c("black","red"),main="In-sample and out-of-sample turbulence")
#' lines(cbind(isturb$epsilon,osturb$epsilon),col=c("black","red"),main="In-sample and out-of-sample epsilon",on=NA)
#' isturb1 <- ISminimax_one(Y, s.k=12)
#' osturb1 <- OSminimax_one(Y, s.k=12, rolling=FALSE, roll.obs=100)
#' plot(cbind(isturb$sturb,isturb1$sturb),col=c("black","red"),main="In-sample turbulence for X and Y")
#' lines(cbind(osturb$sturb,osturb1$sturb),col=c("black","red"),main="Out-of-sample turbulence for X and Y",on=NA,ylim=c(-0.05,1))
#'
#' @import xts
#' @import mice
#' @import fAssets
#' @importFrom zoo rollapplyr
#' @importFrom stats cov
#' @importFrom timeSeries colStats
#' @importFrom stats pf qf
#'
#' @export
OSminimax_one <- function(X, weights=NULL, method=NULL, s.k=1, imp=FALSE, rf=TRUE,
                           rolling=FALSE, roll.obs=100, NOUT=NULL){
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
  # impute missing values
  if (imp){
    x.cn <- colnames(x); colnames(x) <- NULL
    x <- reclass(complete(mice(x,m = 5, method="pmm", maxit=50, seed=500,printFlag = FALSE)),x)
    colnames(x) <- x.cn
  }
  # initialise output vectors
  x <- if (is.vector(X)) matrix(X, ncol=length(X)) else (as.matrix(X))
  m <- dim(x)[1]; n<-dim(x)[2]
  # if NOUT is not specified
  if (is.null(NOUT)){NOUT <- n}
  #
  x.sturb <- x.epsilon <- x[,1]*NA
  # dimension warning
  if((roll.obs-1)<=n){stop("You need more observations to start estimation with (cross-section>number of obs)")}
  #
  if (rf){
    if (s.k==1){ # here we do one future observation
      for (i in roll.obs:m){ #for i>roll.obs we get finally out.of.sample
        if (rolling==FALSE) {x.sel <- x[1:(i-1),]; T <- i} else {x.sel <- x[max(1,i-roll.obs):(i-1),]; T <- roll.obs}
        cols <- which(apply(x.sel,2,function(y) sum(!is.na(y)))>=min(n,roll.obs)) # if rolling window is smaller than no of PFs
        if (length(cols)>0){
          if (!is.null(method)) {
            x.moms <- assetsMeanCov(x.sel[,cols],method=method)
            x.Cov <- getCovRob(x.moms)
            x.mean <- getCenterRob(x.moms)
          }	else {
            x.Cov <- (stats::cov(x.sel[,cols],use="complete.obs")) # even if this deletes observations. But produces correct pos def matrix
            x.mean <- colMeans(x.sel[,cols],na.rm=TRUE)
          }
          x.iCov <- .mpinv(x.Cov)
          x.msw <- colStats(x[max(1,(i-s.k+1)):i,cols,drop=FALSE],FUN=mean,align="right",na.pad=TRUE)-x.mean
          # if there are missing row-values, set them to 0 to be able to matrix multiply below
          x.msw[is.na(x[i,cols])] <- 0
          x.SR <- x.mean %*% x.iCov %*% x.mean
          # NOW NEW: CORRECT TRANSFORMATION
          # x.epsilon[i] <- (1/(T+1) * rowSums((x.msw %*% x.iCov) * x.msw))
          x.epsilon[i] <- NOUT*(T-1)/T/(T-NOUT) * stats::qf(p=stats::pf(q = T*(T-n)/(T+1)/(T-1) * 1/n * rowSums((x.msw %*% x.iCov) * x.msw), df1 = n, df2 = T-n), df1 = NOUT, df2 = T-NOUT)
          x.sturb[i] <- apply(as.matrix(1-(x.epsilon[i]/x.SR)^0.5),1,FUN=function(x){max(x,0)})
        } else {next}
      }
    } else { # for s.k >1 we do subsample (full)
      for (i in roll.obs:m){ #for i>roll.obs we get finally out.of.sample
        if (rolling==FALSE) {x.sel <- x[1:(i),]; T <- i} else {x.sel <- x[max(1,i-roll.obs):(i),]; T <- roll.obs}
        cols <- which(apply(x.sel,2,function(y) sum(!is.na(y)))>=min(n,roll.obs)) # if rolling window is smaller than no of PFs
        if (length(cols)>0 & (T-T/s.k-n+1)>0){
          if (!is.null(method)) {
            x.moms <- assetsMeanCov(x.sel[,cols],method=method)
            x.Cov <- getCovRob(x.moms)
            x.mean <- getCenterRob(x.moms)
          }	else {
            x.Cov <- (stats::cov(x.sel[,cols],use="complete.obs")) # even if this deletes observations. But produces correct pos def matrix
            x.mean <- colMeans(x.sel[,cols],na.rm=TRUE)
          }
          x.iCov <- .mpinv(x.Cov)
          x.msw <- colStats(x[max(1,(i-s.k+1)):i,cols,drop=FALSE],FUN=mean,align="right",na.pad=TRUE)-x.mean
          # if there are missing row-values, set them to 0 to be able to matrix multiply below
          x.msw[is.na(x[i,cols])] <- 0
          x.SR <- x.mean %*% x.iCov %*% x.mean
          # NOW NEW: CORRECT TRANSFORMATION
          # x.epsilon[i] <- s.k*(T-T/s.k-n+1)/(T/s.k-1)/(s.k-1) * (T-1)/T/(T-NOUT) * rowSums((x.msw %*% x.iCov) * x.msw)
          x.epsilon[i] <- NOUT*(T-1)/T/(T-NOUT) * stats::qf(p=stats::pf(q = s.k*(T-T/s.k-n+1)/(T/s.k-1)/(s.k-1) * 1/n * rowSums((x.msw %*% x.iCov) * x.msw), df1 = n, df2 = T-T/s.k-n+1), df1 = NOUT, df2 = T-NOUT)
          #
          x.sturb[i] <- apply(as.matrix(1-(x.epsilon[i]/x.SR)^0.5),1,FUN=function(x){max(x,0)})
        } else {next}
      }
    }
  } else {
    cat("the version without risk-free rate is not yet implemented")
  }
  x.epsilon <- reclass(x.epsilon,try.xts(X))
  x.sturb <- reclass(x.sturb,try.xts(X))
  return(list(sturb=x.sturb,epsilon=x.epsilon))
}
