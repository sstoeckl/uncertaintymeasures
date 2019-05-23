#' @title Calculate portfolio weights for Garlappi model 1.2.2-1.2.3
#'
#' @description Calculate portfolio weights for Garlappi model 1.2.2-1.2.3: Uncertainty about expected returns
#' estimated jointly for all assets and subgroups thereof.
#'
#'
#' @param mu, list of vectors of mus for portfolio subsets
#' @param Sigma, overall covariance matrix
#' @param gamma optional, coefficient of risk aversion
#' @param epsilon, list of coefficients of parameter uncertainty for each subset of parameters defined by mu
#' @param rf optional, should the case with or without the risk-free rate be calculated (only for one PU parameter)
#'
#' @return output of function:
#'     - vector of optimal weights given PU
#'     - sturb: by how much is investment in the optimal portfolio reduced (epsilon vs Sharpe ratio of the other portfolio)
#'
#' @examples
#' require(xts)
#' data(X)
#' mu_in <- list(mu1=colMeans(X[,1:4]),mu2=colMeans(X[,5:7]))
#' Sigma_in <- cov(X[,1:7])
#' gamma <- 3
#' ## First calculate the tangency portfolio for the first 4 assets without parameter uncertainty
#' round(c(1/gamma*.mpinv(Sigma_in[1:4,1:4])%*%mu_in[[1]]),4)
#' # Then do it using the function
#' round(glPF(mu=mu_in[[1]],Sigma=Sigma_in[1:4,1:4],gamma=gamma,epsilon=list(eps1=0),rf=TRUE)$w,4)
#' # Now, we allow for parameter uncertainty for all 4 assets together
#' round(glPF(mu=mu_in[[1]],Sigma=Sigma_in[1:4,1:4],gamma=gamma,epsilon=list(eps1=0.1),rf=TRUE)$w,4)
#' round(glPF(mu=mu_in[[1]],Sigma=Sigma_in[1:4,1:4],gamma=gamma,epsilon=list(eps1=0.2),rf=TRUE)$w,4)
#' # Now we calculate the GMVP
#' round(c(.mpinv(Sigma_in[1:4,1:4])%*%matrix(1,nrow=4,ncol=1)/drop(matrix(1,nrow=1,ncol=4)%*%.mpinv(Sigma_in[1:4,1:4])%*%matrix(1,nrow=4,ncol=1))),4)
#' # and observe how the assets shift towards it
#' round(glPF(mu=mu_in[[1]],Sigma=Sigma_in[1:4,1:4],gamma=gamma,epsilon=list(eps1=0.1),rf=FALSE)$w,4)
#' round(glPF(mu=mu_in[[1]],Sigma=Sigma_in[1:4,1:4],gamma=gamma,epsilon=list(eps1=0.8),rf=FALSE)$w,4)
#' round(glPF(mu=mu_in[[1]],Sigma=Sigma_in[1:4,1:4],gamma=gamma,epsilon=list(eps1=10000),rf=FALSE)$w,4)
#' ## Second: We use two groups for the uncertainty example now
#' # Tangency portfolio for all 7 assets
#' round(c(1/gamma*.mpinv(Sigma_in)%*%unlist(mu_in)),4)
#' # The same using the function
#' round(glPF(mu=mu_in,Sigma=Sigma_in,gamma=gamma,epsilon=list(eps1=0,eps2=0),rf=TRUE)$w,4)
#' # Introduce parameter uncertainty for single assets and groups
#' round(glPF(mu=mu_in,Sigma=Sigma_in,gamma=gamma,epsilon=list(eps1=0.1,eps2=0),rf=TRUE)$w,4)
#' round(glPF(mu=mu_in,Sigma=Sigma_in,gamma=gamma,epsilon=list(eps1=0,eps2=0.1),rf=TRUE)$w,4)
#' round(glPF(mu=mu_in,Sigma=Sigma_in,gamma=gamma,epsilon=list(eps1=0.1,eps2=0.05),rf=TRUE)$w,4)
#' ## Third: We do it for single assets (two)
#' mu_in <- list(mu1=colMeans(X[,1]),mu2=colMeans(X[,2]))
#' Sigma_in <- cov(X[,1:2])
#' gamma <- 3
#' # Tangency portfolio for 2 assets
#' round(c(1/gamma*.mpinv(Sigma_in)%*%unlist(mu_in)),4)
#' # The same using the function
#' round(glPF(mu=mu_in,Sigma=Sigma_in,gamma=gamma,epsilon=list(eps1=0,eps2=0),rf=TRUE)$w,4)
#' # Introduce parameter uncertainty for single assets and groups
#' round(glPF(mu=mu_in,Sigma=Sigma_in,gamma=gamma,epsilon=list(eps1=0.1,eps2=0),rf=TRUE)$w,4)
#' round(glPF(mu=mu_in,Sigma=Sigma_in,gamma=gamma,epsilon=list(eps1=0,eps2=0.1),rf=TRUE)$w,4)
#' round(glPF(mu=mu_in,Sigma=Sigma_in,gamma=gamma,epsilon=list(eps1=0.0002,eps2=0.01),rf=TRUE)$w,4)
#' # Sturb
#' round(glPF(mu=mu_in,Sigma=Sigma_in,gamma=gamma,epsilon=list(eps1=0.1,eps2=0),rf=TRUE)$sturb,4)
#' round(glPF(mu=mu_in,Sigma=Sigma_in,gamma=gamma,epsilon=list(eps1=0,eps2=0.1),rf=TRUE)$sturb,4)
#' round(glPF(mu=mu_in,Sigma=Sigma_in,gamma=gamma,epsilon=list(eps1=0.0002,eps2=0.01),rf=TRUE)$sturb,4)
#'
#' @import nleqslv
#' @import xts
#'
#' @keywords internal
#'
#' @export
glPF <- function(mu, Sigma, gamma=3, epsilon, rf){
  if (is.list(mu)){n <- length(mu)} else {mu <- list(mu); n<-1}
  N <- length(unlist(mu))
  if (n==1){
    if(rf){
      # with risk-free rate
      w <- c(max(1-epsilon[[1]]/drop(mu[[1]]%*%.mpinv(Sigma)%*%mu[[1]]),0)*1/gamma*.mpinv(Sigma)%*%mu[[1]])
    } else {
      # without risk-free rate
      A <- drop(rep(1,ncol(Sigma))%*%.mpinv(Sigma)%*%rep(1,ncol(Sigma)))
      B <- drop(mu[[1]]%*%.mpinv(Sigma)%*%rep(1,ncol(Sigma)))
      C <- drop(mu[[1]]%*%.mpinv(Sigma)%*%mu[[1]])
      d_sigma_p <- function(sigma_p,A,B,C,gamma,epsilon){
        return(A*gamma^2*sigma_p^4 +
                 2*A*gamma*sqrt(epsilon[[1]])*sigma_p^3 +
                 (A*epsilon[[1]] - A*C + B^2-gamma^2)*sigma_p^2 - 2*gamma*sqrt(epsilon[[1]])*sigma_p-epsilon[[1]]-sigma_p)
      }
      sigma_p <- nleqslv::nleqslv(0.1,d_sigma_p,A=A,B=B,C=C,gamma=gamma,epsilon=epsilon[[1]])$x
      # weights
      w <- c(sigma_p/(sqrt(epsilon[[1]]) + gamma*sigma_p)*.mpinv(Sigma)%*%t(mu[[1]]-1/A*(B-(sqrt(epsilon[[1]])+gamma*sigma_p)/sigma_p)*matrix(1,ncol=length(mu[[1]]),nrow=1)))
    }
  } else {
    if(rf){
      #cat("Calculating PU optimal portfolio for",n,"groups and",N,"assets.\n")
      # define helper function
      sub_fn <- function(mu_in,mu_out,Sigma_in,Sigma_out,eps_in,w_out,gamma){
        g <- function(w_out,mu_in,gamma,Sigma_out){
          return(mu_in - gamma * (Sigma_out %*% w_out))
        }
        out <- 1/gamma * .mpinv(Sigma_in) %*% g(w_out,mu_in,gamma,Sigma_out) *
          max(1-sqrt(eps_in/t(g(w_out,mu_in,gamma,Sigma_out)) %*% .mpinv(Sigma_in) %*% g(w_out,mu_in,gamma,Sigma_out)),0)
        return(out)
      }
      # now create optimizer function
      fn <- function(w,mu,Sigma,epsilon,gamma){
        n <- length(mu); nvec <- plyr::laply(mu,length)
        N <- length(unlist(mu))
        dw <- NULL
        for (i in 1:n){
          vec <- seq(from=max(nvec[i-1],0)+1,to=sum(nvec[1:i]))
          temp_dw <- (c(sub_fn(mu_in=mu[[i]],
                               mu_out=unlist(mu[[-i]]),
                               Sigma_in=Sigma[vec,vec],
                               Sigma_out=Sigma[vec,-vec],
                               w_out=w[-vec],
                               eps_in=epsilon[[i]],gamma=gamma)) - w[vec])
          dw <- c(dw,temp_dw)
        }
        return(dw)
      }
      #fn(rep(1/N,N),mu,Sigma,epsilon,gamma)
      # now optimize
      w <- nleqslv::nleqslv(rep(1/N,N), fn,
                            mu=mu,
                            Sigma=Sigma,
                            epsilon=epsilon,
                            gamma=gamma)$x
    } else {stop("version without risk-free rate is not implemented yet!")}
  }
  # sturb: Similar to before, see how far this moves the asset out of the market in relation to its former weights
  sub_sturb <- function(mu_in,mu_out,Sigma_in,Sigma_out,eps_in,w_out,gamma){
    g <- function(w_out,mu_in,gamma,Sigma_out){
      return(mu_in - gamma * (Sigma_out %*% w_out))
    }
    out <- max(1-sqrt(eps_in/t(g(w_out,mu_in,gamma,Sigma_out)) %*% .mpinv(Sigma_in) %*% g(w_out,mu_in,gamma,Sigma_out)),0)
    return(out)
  }
  n <- length(mu); nvec <- plyr::laply(mu,length)
  N <- length(unlist(mu))
  sturb <- NULL
  for (i in 1:n){
    vec <- seq(from=max(nvec[i-1],0)+1,to=sum(nvec[1:i]))
    temp_sturb <- (c(sub_sturb(mu_in=mu[[i]],
                               mu_out=unlist(mu[[-i]]),
                               Sigma_in=Sigma[vec,vec],
                               Sigma_out=Sigma[vec,-vec],
                               w_out=w[-vec],
                               eps_in=epsilon[[i]],gamma=gamma)))
    sturb <- c(sturb,temp_sturb)
  }
  names(w) <- names(unlist(mu))
  names(sturb) <- names(unlist(factor))
  return(list(w=w,sturb=sturb))
}

#' @title Calculate portfolio weights for Garlappi model 1.2.3
#'
#' @description Calculate portfolio weights for Garlappi model 1.2.2-1.2.3: Uncertainty about expected returns
#' estimated jointly for all assets and subgroups thereof.
#'
#' @param factors, list of factor premia xts
#' @param PFs, list of portfolios (xts) to calculate PU for each factor premium
#' @param gamma optional, coefficient of risk aversion (standard = 3)
#' @param s.k optional, lookback window for current observations that are related to the long-term mean
#' @param imp optional, should missing values be imputed? (standard: FALSE)
#' @param rolling.PF optional, should portfolio moments be calculated based on a rolling window?
#' @param rolling.eps optional, should uncertainty moments be calculated based on a rolling window?
#' @param roll.obs.PF optional, rolling.PF=TRUE: Length of rolling window; rolling.PF=FALSE: Initial estimation window size
#' @param roll.obs.eps optional, rolling.eps=TRUE: Length of rolling window; rolling.eps=FALSE: Initial estimation window size
#' @param rf optional, should the case with or without the risk-free rate be calculated (only for one PU parameter)
#'
#' @return
#'
#' @examples
#' data(factors)
#' data(X)
#' data(Y)
#' # out-of-sample exercise
#' factors <- list(fact1=factors[,1],fact2=factors[,2])
#' PFs <- list(PF1=X,PF2=Y)
#' w1 <- OSglPF(factors,PFs, s.k=12, imp=FALSE, rolling.PF=FALSE, rolling.eps=FALSE, roll.obs.PF=120, roll.obs.eps=120)
#' plot(w1$epsilon)
#' plot(w1$weights)
#' plot(w1$sturb)
#' w2 <- OSglPF(factors,PFs, s.k=12, imp=FALSE, rolling.PF=TRUE, rolling.eps=TRUE, roll.obs.PF=120, roll.obs.eps=120)
#' plot(w2$epsilon)
#' plot(w2$weights)
#' plot(w2$sturb)
#' w3 <- OSglPF(factors,PFs, s.k=12, imp=FALSE, rolling.PF=FALSE, rolling.eps=TRUE, roll.obs.PF=120, roll.obs.eps=120)
#' plot(w3$epsilon)
#' plot(w3$weights)
#' plot(w3$sturb)
#'
#' @importFrom plyr ldply
#' @import xts
#' @import zoo
#'
#' @keywords internal
#'
#' @export
OSglPF <- function(factors,PFs, s.k=1, imp=FALSE, rolling.PF=FALSE, rolling.eps=FALSE, roll.obs.PF=120,roll.obs.eps=120, rf=TRUE, gamma=3){
  # method = c("cov", "mve", "mcd", "MCD", "OGK", "nnve", "shrink", "bagged")
  if (!requireNamespace("xts", quietly = TRUE)) {
    stop("Package \"xts\" needed for this function to work. Please install it.",call. = FALSE)}
  if (!requireNamespace("fAssets", quietly = TRUE)) {
    stop("Package \"fAssets\" needed for this function to work. Please install it.",call. = FALSE)}
  if(imp){
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("Package \"mice\" needed for this function to work. Please install it.",call. = FALSE)}
  }
  ## Function Main
  N <- length(factors)
  epsilon <- list()
  # 0. Adjust factors and all datasets if necessary
  min_date1 <- max(plyr::ldply(PFs,function(x){as.Date(time(x)[1])})$V1)
  min_date2 <- max(plyr::ldply(factors,function(x){as.Date(time(x)[1])})$V1)
  PFs <- plyr::llply(PFs,function(x){x[paste0(min_date1,"/")]})
  # 1. Now calculate OS sturb and save epsilons
  for (i in 1:N){
    epsilon[[i]] <- OSminimax_one(X = PFs[[i]], s.k = s.k, rolling=rolling.eps, roll.obs = roll.obs.eps, NOUT=N)$epsilon
  }
  min_date_eps <- min(plyr::ldply(epsilon,function(x){as.Date(time(na.omit(x))[1])})$V1)
  stopifnot(as.yearmon(min_date_eps)<=(as.yearmon(min_date2)+roll.obs.PF/12))
  # 2. Now loop through each point in time and create portfolio weights
  T <- nrow(factors[[1]])
  weights <- matrix(NA,nrow=T,ncol=N)
  sturb <- matrix(NA,nrow=T,ncol=N)
  for (j in roll.obs.PF:T){
    if (rolling.PF==FALSE) {
      mu_in <- plyr::llply(factors,function(x,j){mean(x[1:(j-1)])},j=j)
      Sigma_in <- cov(t(plyr::laply(factors,function(x,j){x[1:(j-1)]},j=j)))
      epsilon_in <- plyr::llply(epsilon,function(x,j){x[j]},j=j)
      temp <- glPF(mu=mu_in, Sigma=Sigma_in, gamma=gamma, epsilon=epsilon_in, rf=rf)
      weights[j,] <- temp$w
      sturb[j,] <- temp$sturb
    } else {
      mu_in <- plyr::llply(factors,function(x,j,roll.obs){mean(x[max(1,j-roll.obs):(j-1)])},j=j,roll.obs=roll.obs.PF)
      Sigma_in <- cov(t(plyr::laply(factors,function(x,j,roll.obs){x[max(1,j-roll.obs):(j-1)]},j=j,roll.obs=roll.obs.PF)))
      epsilon_in <- plyr::llply(epsilon,function(x,j,roll.obs){mean(x[j])},j=j)
      temp <- glPF(mu=mu_in, Sigma=Sigma_in, gamma=gamma, epsilon=epsilon_in, rf=rf)
      weights[j,] <- temp$w
      sturb[j,] <- temp$sturb
    }
  }
  weights <- reclass(weights,try.xts(factors[[1]]))
  sturb <- reclass(sturb,try.xts(factors[[1]]))
  epsilon.out <- Reduce(cbind,epsilon)
  return(list(weights=weights,epsilon=epsilon.out,sturb=sturb))
}

