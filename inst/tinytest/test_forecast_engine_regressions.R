# Regime forecasts must follow connected Markov paths.
set.seed(20260806)
S             = 10000L
transition    = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE)
sigma2        = array(rep(c(1, 2), S), c(1, 2, S))
transitions   = array(rep(transition, S), c(2, 2, S))
terminal_state = matrix(rep(c(1, 0), S), 2, S)

msh_paths = .Call(
  "_bsvars_forecast_sigma2_msh",
  sigma2,
  transitions,
  terminal_state,
  2L,
  PACKAGE = "bsvars"
)

expect_equal(
  mean(msh_paths[1, 1, ] == 1 & msh_paths[1, 2, ] == 1),
  0.81,
  tolerance = 0.02,
  info = "forecast_sigma2_msh: samples a connected Markov path."
)
expect_true(
  cor(msh_paths[1, 1, ], msh_paths[1, 2, ]) > 0.5,
  info = "forecast_sigma2_msh: preserves serial state dependence."
)

set.seed(20260806)
hmsh_transitions = replicate(
  S,
  array(transition, c(2, 2, 1)),
  simplify = FALSE
)
hmsh_terminal_state = array(rep(c(1, 0), S), c(2, 1, S))

hmsh_paths = .Call(
  "_bsvars_forecast_sigma2_hmsh",
  sigma2,
  hmsh_transitions,
  hmsh_terminal_state,
  2L,
  PACKAGE = "bsvars"
)

expect_equal(
  mean(hmsh_paths[1, 1, ] == 1 & hmsh_paths[1, 2, ] == 1),
  0.81,
  tolerance = 0.02,
  info = "forecast_sigma2_hmsh: samples a connected Markov path."
)
expect_true(
  cor(hmsh_paths[1, 1, ], hmsh_paths[1, 2, ]) > 0.5,
  info = "forecast_sigma2_hmsh: preserves serial state dependence."
)


# Gaussian sampling failures must stop the forecast.
expect_error(
  .Call(
    "_bsvars_forecast_bsvars",
    array(1, c(1, 1, 1)),
    array(0, c(1, 1, 1)),
    array(-1, c(1, 1, 1)),
    0,
    matrix(NA_real_, 1, 1),
    matrix(NA_real_, 1, 1),
    1L,
    PACKAGE = "bsvars"
  ),
  info = "forecast_bsvars: propagates covariance and sampling failures."
)
expect_error(
  .Call(
    "_bsvars_forecast_bsvars",
    array(diag(2), c(2, 2, 1)),
    array(0, c(2, 2, 1)),
    array(c(-1, 1), c(2, 1, 1)),
    c(0, 0),
    matrix(NA_real_, 1, 1),
    matrix(c(NA_real_, 0), 1, 2),
    1L,
    PACKAGE = "bsvars"
  ),
  info = "forecast_bsvars: propagates conditional sampling failures."
)


# Conditional forecasts must return conditional moments.
unconditional_covariance = matrix(c(1, 0.5, 0.5, 1), 2, 2)
B_inverse                = t(chol(unconditional_covariance))

conditional_draw = .Call(
  "_bsvars_mvnrnd_cond",
  c(2, NA_real_),
  c(0, 0),
  unconditional_covariance,
  PACKAGE = "bsvars"
)
expect_equal(
  conditional_draw[1],
  2,
  info = "mvnrnd_cond: retains its vector-returning public interface."
)

conditional_output       = .Call(
  "_bsvars_forecast_bsvars",
  array(solve(B_inverse), c(2, 2, 1)),
  array(0, c(2, 2, 1)),
  array(1, c(2, 1, 1)),
  c(0, 0),
  matrix(NA_real_, 1, 1),
  matrix(c(2, NA_real_), 1, 2),
  1L,
  PACKAGE = "bsvars"
)

expect_equal(
  conditional_output$forecasts[1, 1, 1],
  2,
  info = "forecast_bsvars: respects the fixed coordinate."
)
expect_equal(
  conditional_output$forecast_mean[, 1, 1],
  c(2, 1),
  tolerance = 1e-10,
  info = "forecast_bsvars: returns the full conditional mean."
)
expect_equal(
  conditional_output$forecast_cov[[1]][, , 1],
  matrix(c(0, 0, 0, 0.75), 2, 2),
  tolerance = 1e-10,
  info = "forecast_bsvars: returns the full conditional covariance."
)
