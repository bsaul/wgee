#' Make cluster level matrices necessary for computations
#'
#' @param df the data.frame in \code{\link{extract_geeglm}}. It must contain
#' the following variables in the following order: \code{response}, \code{fitted}, \code{weights}.
#' Remaining (covariate) variables must begin at column index 4.
#'
#' @export

make_cluster_matrices <- function(df){

  #TODO: DI and Vi.inv are functions of the family and link function
  # These are currently set up for Gaussian family with identity link function ONLY!
  list(ri = as.matrix(df$response - df$fitted),
       Di = as.matrix(df[, 4:ncol(df)]),
       Wi = diag(df$weights),
       Vi.inv = diag(1, nrow = nrow(df), ncol = nrow(df)) )
}

#' Cluster-level O matrix
#'
#' Computes
#' \deqn{ D_i^T V_i^{-1} W_i D_i}
#'
#' @param D matrix of derivatives of mean function
#' @param W weight matrix
#' @param V.inv inverse variance matrix
#'
#' @export

Oi_matrix <- function(D, W, V.inv){
  t(D) %*% V.inv %*% W %*% D
}

#' Sums a list of matrices
#'
#' @param matrix_list list of matrices to sum across (must all be of same dimensions)
#' @param geeglm_ingredients result of \code{\link{extract_geeglm}}
#' @importFrom magrittr "%>%"
#'

summerizer <- function(matrix_list, geeglm_ingredients, solve = FALSE){
  p <- geeglm_ingredients$ncoef
  N <- geeglm_ingredients$nclust

  out <-  matrix_list %>%
    unlist() %>%
    array(dim = c(p, p, N)) %>%
    apply(1:2, sum)

  if(solve == TRUE){
    solve(out)
  } else {
    out
  }
}

#' Bread maker
#'
#' Computes
#' \deqn{ \sum_{i = 1}^N D_i^T V_i^{-1} W_i D_i}
#'
#' @param matrix_list list of matrices
#' @param geeglm_ingredients result of \code{\link{extract_geeglm}}
#'
#' @export

bread <- function(matrix_list, geeglm_ingredients){

  lapply(matrix_list, function(m) wgee::Oi_matrix(D = m$Di, W = m$Wi, V.inv = m$Vi.inv ) ) %>%
    summerizer(geeglm_ingredients, solve = TRUE)
}


#' Cluster-level U matrix
#'
#' Computes
#' \deqn{ D_i^T V_^{-1} W_i r_i}
#'
#' @param D matrix of derivatives of mean function
#' @param V.inv inverse variance matrix
#' @param W matrix of weight
#' @param r matrix of residuals
#'
#' @export

Ui_matrix <- function(D, V.inv, W, r){
  t(D) %*% V.inv %*% W %*% r
}


#' Cluster-level H matrix
#'
#' Computes
#' \deqn{ D_i O D_i^T V_^{-1}}
#'
#' where \eqn{O} is the result of the \code{\link{bread}} function
#'
#' @param D matrix of derivatives of mean function
#' @param V.inv inverse variance matrix
#' @param bread result of the \code{\link{bread}} function
#'
#' @export

Hi_matrix <- function(D, V.inv, bread){
  D %*% bread %*% t(D) %*% V.inv
}

#' Appends a list of Ui and Hi matrices
#'
#' @param matrix_list
#' @param bread result of the \code{\link{bread}} function
#'
#' @export

add_U_H <- function(matrix_list, bread){
  lapply(matrix_list, FUN = function(m) {
    append(m, list(
      Ui = wgee::Ui_matrix(D = m$Di, V.inv = m$Vi.inv, W = m$Wi, r = m$ri),
      Hi = wgee::Hi_matrix(D = m$Di, V.inv = m$Vi.inv, bread) ) )
  } )
}

#' Meat maker
#'
#' Computes
#' \deqn{ \sum_{i = 1}^N U_i U_i^T}
#'
#' @param matrix_list list of matrices
#' @param geeglm_ingredients result of \code{\link{extract_geeglm}}
#'
#' @export


meat <- function(matrix_list, geeglm_ingredients){
  lapply(matrix_list, function(m) m$Ui %*% t(m$Ui)) %>%
    summerizer(geeglm_ingredients, solve = FALSE)
}
