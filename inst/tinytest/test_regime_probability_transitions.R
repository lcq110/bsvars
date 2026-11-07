transition = matrix(c(0.8, 0.2, 0.1, 0.9), 2, 2, byrow = TRUE)
initial_probability = c(0.6, 0.4)

common_posterior = list(
  B = array(1, c(1, 1, 1)),
  A = array(0, c(1, 1, 1)),
  sigma2 = array(c(1, 1), c(1, 2, 1)),
  PR_TR = array(transition, c(2, 2, 1)),
  pi_0 = matrix(initial_probability, 2, 1)
)
Y = matrix(0, 1, 1)
X = matrix(0, 1, 1)

common_filtered = .Call(
  "_bsvars_bsvars_filter_forecast_smooth",
  common_posterior,
  Y,
  X,
  FALSE,
  FALSE,
  PACKAGE = "bsvars"
)
common_forecasted = .Call(
  "_bsvars_bsvars_filter_forecast_smooth",
  common_posterior,
  Y,
  X,
  TRUE,
  FALSE,
  PACKAGE = "bsvars"
)

expect_equal(
  as.numeric(common_forecasted[, 1, 1]),
  as.numeric(t(transition) %*% common_filtered[, 1, 1]),
  tolerance = 1e-12,
  info = "Common-state forecast probabilities sum incoming transition probabilities."
)
expect_equal(
  sum(common_forecasted[, 1, 1]),
  1,
  tolerance = 1e-12,
  info = "Common-state forecast probabilities remain normalized."
)

N = 2
T = 2
hmsh_posterior = list(
  B = array(diag(N), c(N, N, 1)),
  A = array(0, c(N, 1, 1)),
  sigma2 = array(matrix(c(1, 4, 1, 4), N, 2, byrow = TRUE), c(N, 2, 1)),
  PR_TR_cpp = list(array(rep(transition, N), c(2, 2, N))),
  pi_0 = array(rep(initial_probability, N), c(2, N, 1))
)
X = matrix(0, 1, T)
Y_1 = matrix(c(1, 1, 0, 0), N, T, byrow = TRUE)
Y_2 = matrix(c(1, 1, 20, -20), N, T, byrow = TRUE)

hmsh_filtered_1 = .Call(
  "_bsvars_bsvars_filter_forecast_smooth_hmsh",
  hmsh_posterior,
  Y_1,
  X,
  FALSE,
  FALSE,
  PACKAGE = "bsvars"
)[[1]]
hmsh_filtered_2 = .Call(
  "_bsvars_bsvars_filter_forecast_smooth_hmsh",
  hmsh_posterior,
  Y_2,
  X,
  FALSE,
  FALSE,
  PACKAGE = "bsvars"
)[[1]]

expect_equal(
  hmsh_filtered_1[, , 1],
  hmsh_filtered_2[, , 1],
  tolerance = 1e-12,
  info = "Shock-specific HMSH probabilities depend only on that shock's observations and variances."
)

hmsh_forecasted = .Call(
  "_bsvars_bsvars_filter_forecast_smooth_hmsh",
  hmsh_posterior,
  Y_1,
  X,
  TRUE,
  FALSE,
  PACKAGE = "bsvars"
)[[1]]

expect_equal(
  hmsh_forecasted[, , 1],
  t(transition) %*% hmsh_filtered_1[, , 1],
  tolerance = 1e-12,
  info = "HMSH forecast probabilities sum incoming transition probabilities shock by shock."
)
expect_equal(
  colSums(hmsh_forecasted[, , 1]),
  rep(1, T),
  tolerance = 1e-12,
  info = "HMSH forecast probabilities remain normalized."
)
