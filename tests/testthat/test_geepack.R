library(geepack)
library(geesmv)

data("ohio")


#-----------------------------------------------------------------------------------#
# Gaussian with identity link
#-----------------------------------------------------------------------------------#

context("Gaussian, identity link")

geefit <- geepack::geeglm(resp ~ age, family = gaussian, data = ohio, id = id)
wgeefit_none <- wgee::geeglm_w(resp ~ age, family = gaussian, data = ohio, id = id,
                               correction = 'none')

geesmv_md <- geesmv::GEE.var.md(resp ~ age, family = gaussian, data = ohio, id = id)
wgeefit_md <- wgee::geeglm_w(resp ~ age, family = gaussian, data = ohio, id = id,
                             correction = 'MD')

geesmv_kc <- geesmv::GEE.var.kc(resp ~ age, family = gaussian, data = ohio, id = id)
wgeefit_kc <- wgee::geeglm_w(resp ~ age, family = gaussian, data = ohio, id = id,
                            correction = 'KC')

test_that("wgee returns same coefficients and standard errors as geepack when no correction used", {
  expect_equal(summary(geefit)$coefficients, summary(wgeefit_none)$coefficients)
})


test_that("wgee returns same covariance matrix as geesmv when MD correction used ", {
  test1 <- geesmv_md$cov.beta
  names(test1) <- NULL
  test2 <- diag(wgeefit_md$geese$vbeta)
  expect_equal(test1, test2)
})

test_that("wgee returns same covariance matrix as geesmv when KC correction used ", {
  test1 <- geesmv_kc$cov.beta
  names(test1) <- NULL
  test2 <- diag(wgeefit_kc$geese$vbeta)
  expect_equal(test1, test2)
})
