#' @title Helper function for pseudoinverse
#'
#' @description Calculate the Penrose-Moore pseudoinverse based on adjusting the eigenvalues of the matrix
#'
#' @param A, Any matrix that you wish to invert
#' @param eps optional, the level of tolerance to use
#'
#' @examples
#' A1 <- matrix(2,1,1)
#' .mpinv(A1)
#' A2 <- 3
#' .mpinv(A2)
#' A3 <- matrix(c(2,0.1,0.3,0.1,1,0.5),2,3)
#' .mpinv(A3)
#' round(A3%*%.mpinv(A3),2)
#'
#' @export
.mpinv <- function(A, eps = 1e-13) {
  if(!is.matrix(A)){
    return(matrix(1/A,1,1))
  } else if(min(dim(A))==1) {
    return(matrix(1/A,1,1))
  } else {
    s <- svd(A)
    e <- s$d
    e[e > eps] <- 1/e[e > eps]
    return(s$v %*% diag(e) %*% t(s$u))
  }
}
