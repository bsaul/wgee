#' geeglm_w
#'
#' @param geeobj
#' @param Ufun
#'
#' @export
#'

geeglm_w <- function(formula, data, id, family = gaussian,
                     correction = 'none', ...)
{
  geecall <- match.call()
  geecall[[1]] <- quote(geepack::geeglm)
  geecall[['correction']] <- NULL

  fit <- eval(geecall)


  if (correction == 'MD'){
    Ufun = wgee::U_MD
  } else if (correction == 'KC'){
    Ufun = wgee::U_KC
  } else {
    Ufun = wgee::U_standard
  }

  # Switch out covariance matrices
  fit$geese$vbeta <- wgee::sandwich(fit, Ufun = Ufun)

  return(fit)
}
