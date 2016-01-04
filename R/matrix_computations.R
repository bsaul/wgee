#' Taken from: http://www4.stat.ncsu.edu/~li/software/SparseSDR.R
#' square-root-inverse of a matrix
#' @export

sqrt_matrix <-function(M)
{
  ei <-eigen(M)
  d  <-ei$values
  d  <-(d+abs(d))/2
  d2 <-1 / sqrt(d)
  d2[d < 0] <- 0
  ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(ans)
}
