transition = matrix(c(0.8, 0.2, 0.1, 0.9), 2, 2, byrow = TRUE)
initial_probability = c(0.6, 0.4)
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
