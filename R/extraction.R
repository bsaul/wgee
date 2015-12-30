#' Extract necessary elements from geeglm object
#'
#' @param geeobj
#'
#' @export

extract_geeglm <- function(geeobj){
  if(!('geeglm' %in% class(geeobj) ) ){
    stop("Object must be of class geeglm")
  }

  out <- list(family = geeobj$family$family,
              link   = geeobj$family$link,
              beta   = coef(geeobj),
              ncoef  = length(coef(geeobj)),
              nclust = length(unique(geeobj$id)),
              id     = geeobj$id,
              data   = cbind(data.frame(response = geeobj$y,
                                        fitted   = fitted(geeobj),
                                        weights  = geeobj$weights),
                             model.matrix(geeobj) ) )
 return(out)
}
