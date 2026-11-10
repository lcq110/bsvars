irf = array(0, c(2, 2, 2))
irf[, , 1] = matrix(c(1, 1, 0, 1), 2, 2, byrow = TRUE)
irf[, , 2] = matrix(c(2, 1, 0, 1), 2, 2, byrow = TRUE)
forecast_sigma2 = array(c(4, 1), c(2, 1, 1))
sigma2_T = matrix(c(1, 1), 2, 1)

fevd = .Call(
  "_bsvars_bsvars_fevd_heterosk",
  list(irf),
  forecast_sigma2,
  sigma2_T,
  PACKAGE = "bsvars"
)[[1]]

expect_equal(
  as.numeric(fevd[1, , 2]),
  c(80, 20),
  tolerance = 1e-12,
  info = "Heteroskedastic FEVD pairs each future variance with the matching reverse-horizon response."
)
