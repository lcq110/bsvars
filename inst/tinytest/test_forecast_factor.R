call_forecast_engine = function(B, sigma2, horizon = dim(sigma2)[2]) {
  N = dim(B)[1]
  S = dim(B)[3]

  .Call(
    "_bsvars_forecast_bsvars",
    B,
    array(0, c(N, N + 1L, S)),
    sigma2,
    rep(0, N + 1L),
    matrix(NA_real_, horizon, 1L),
    matrix(NA_real_, horizon, N),
    as.integer(horizon),
    PACKAGE = "bsvars"
  )
}

log_density_chol = function(error, covariance) {
  covariance_chol = chol(covariance)
  standardised_error = forwardsolve(t(covariance_chol), error)

  -0.5 * (
    length(error) * log(2 * pi) +
      2 * sum(log(diag(covariance_chol))) +
      sum(standardised_error^2)
  )
}

log_density_factor = function(error, B, sigma2) {
  log_abs_det_B = as.numeric(determinant(B, logarithm = TRUE)$modulus)
  structural_error = drop(B %*% error)

  log_abs_det_B - 0.5 * (
    length(error) * log(2 * pi) +
      sum(log(sigma2)) +
      sum(structural_error^2 / sigma2)
  )
}


# The structural and covariance representations agree when both are stable.
well_B = matrix(
  c(
    1.0, 0.2, -0.1,
    0.1, 1.1, 0.05,
    -0.2, 0.1, 0.9
  ),
  3L,
  byrow = TRUE
)
well_sigma2 = c(0.5, 1.2, 2.0)
well_output = call_forecast_engine(
  array(well_B, c(3L, 3L, 1L)),
  array(well_sigma2, c(3L, 1L, 1L))
)
well_covariance = well_output$forecast_cov[[1L]][, , 1L]
well_error = c(0.3, -0.4, 0.2)

expect_identical(
  names(well_output),
  c("forecasts", "forecast_mean", "forecast_cov"),
  info = "forecast_bsvars keeps its seven-argument, three-field raw interface."
)
expect_equal(
  well_B %*% well_covariance %*% t(well_B),
  diag(well_sigma2),
  tolerance = 1e-12,
  info = "The returned covariance remains consistent with B and sigma2."
)
expect_equal(
  log_density_factor(well_error, well_B, well_sigma2),
  log_density_chol(well_error, well_covariance),
  tolerance = 1e-12,
  info = "Structural-factor and legacy Cholesky Gaussian scores agree."
)


# One standard-normal vector is consumed for every posterior draw and horizon.
rng_N = 2L
rng_H = 3L
rng_S = 2L
rng_B = array(0, c(rng_N, rng_N, rng_S))
rng_B[, , 1L] = diag(c(2, 4))
rng_B[, , 2L] = diag(c(0.5, 5))
rng_sigma2 = array(
  c(
    1, 4,
    4, 1,
    9, 16,
    16, 9,
    4, 25,
    1, 4
  ),
  c(rng_N, rng_H, rng_S)
)

set.seed(20260810)
expected_forecasts = array(NA_real_, c(rng_N, rng_H, rng_S))
for (s in seq_len(rng_S)) {
  for (h in seq_len(rng_H)) {
    one_step = call_forecast_engine(
      array(rng_B[, , s], c(rng_N, rng_N, 1L)),
      array(rng_sigma2[, h, s], c(rng_N, 1L, 1L))
    )
    expected_forecasts[, h, s] = one_step$forecasts[, 1L, 1L]
  }
}
expected_rng_state = .Random.seed

set.seed(20260810)
rng_output = call_forecast_engine(rng_B, rng_sigma2, rng_H)
actual_rng_state = .Random.seed

expect_equal(
  rng_output$forecasts,
  expected_forecasts,
  tolerance = 1e-14,
  info = "Unconditional forecasts use one N-vector of standard normals per s and h."
)
expect_identical(
  actual_rng_state,
  expected_rng_state,
  info = "The joint call consumes the same RNG as one N-vector per s and h."
)
set.seed(20260810)
invisible(rnorm(rng_N * rng_H * rng_S))
expect_identical(
  actual_rng_state,
  .Random.seed,
  info = "Unconditional forecasts consume exactly N RNG values per s and h."
)


# Squaring an ill-conditioned factor can erase its small singular direction.
delta = 2^-30
ill_B = matrix(
  c(0.5, 0.5, 1 / (2 * delta), -1 / (2 * delta)),
  2L,
  byrow = TRUE
)
ill_warnings = character()
set.seed(481516)
ill_reference = call_forecast_engine(
  array(diag(2), c(2L, 2L, 1L)),
  array(1, c(2L, 1L, 1L))
)$forecasts[, 1L, 1L]
set.seed(481516)
ill_output = withCallingHandlers(
  call_forecast_engine(
    array(ill_B, c(2L, 2L, 1L)),
    array(1, c(2L, 1L, 1L))
  ),
  warning = function(warning) {
    ill_warnings <<- c(ill_warnings, conditionMessage(warning))
    invokeRestart("muffleWarning")
  }
)
ill_draw = ill_output$forecasts[, 1L, 1L]
ill_covariance = ill_output$forecast_cov[[1L]][, , 1L]
ill_error = c(0.25, -0.25)

expect_identical(
  ill_warnings,
  character(),
  info = "The factor solve does not warn about an approximate solution."
)
expect_error(
  chol(ill_covariance),
  info = "The pinned ill-conditioned covariance has no usable Cholesky factor."
)
expect_true(
  all(is.finite(ill_draw)),
  info = "The structural-factor draw remains finite when covariance Cholesky fails."
)
expect_equal(
  drop(ill_B %*% ill_draw),
  ill_reference,
  tolerance = 1e-6,
  info = "The ill-conditioned draw satisfies the structural system."
)
expect_true(
  is.finite(log_density_factor(ill_error, ill_B, c(1, 1))),
  info = "The structural Gaussian score remains finite without covariance repair."
)


# At machine precision, no approximate solve is substituted for an exact failure.
near_singular_B = matrix(
  c(1, 1, 1, 1 + 2^-52),
  2L,
  byrow = TRUE
)
near_singular_warnings = character()
near_singular_result = tryCatch(
  withCallingHandlers(
    call_forecast_engine(
      array(near_singular_B, c(2L, 2L, 1L)),
      array(1, c(2L, 1L, 1L))
    ),
    warning = function(warning) {
      near_singular_warnings <<- c(
        near_singular_warnings,
        conditionMessage(warning)
      )
      invokeRestart("muffleWarning")
    }
  ),
  error = identity
)

near_singular_valid = if (inherits(near_singular_result, "error")) {
  nzchar(conditionMessage(near_singular_result))
} else {
  near_singular_draw = near_singular_result$forecasts[, 1L, 1L]
  all(is.finite(near_singular_draw)) &&
    all(is.finite(near_singular_B %*% near_singular_draw))
}

expect_true(
  near_singular_valid,
  info = "A near-singular solve either returns a finite residual or fails explicitly."
)
expect_false(
  any(grepl("attempting approx solution", near_singular_warnings, fixed = TRUE)),
  info = "A near-singular factor solve never attempts an approximate solution."
)
