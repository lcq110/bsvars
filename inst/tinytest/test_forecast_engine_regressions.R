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
