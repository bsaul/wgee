#' Extract necessary elements from geeglm object
#'
#' @param geeobj
#'
#' @export
#

extract_geeglm <- function(geeobj){
  if(typeof(geeobj) != 'geeglm'){
    stop("Object must be of type geeglm")
  }
}
