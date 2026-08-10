
data(us_fiscal_lsuw)

check_forecast_factor = function(forecasts, posterior, horizon, model) {
  posterior_B = posterior$posterior$B
  factor = forecasts$forecast_factor
  expected_dimensions = c(
    dim(posterior_B)[1],
    as.integer(horizon),
    dim(posterior_B)[3]
  )

  expect_identical(
    names(forecasts)[1:5],
    c("forecasts", "forecast_mean", "forecast_cov", "forecast_covariance", "Y"),
    info = paste(model, "keeps the existing first five forecast fields.")
  )
  expect_identical(
    names(factor),
    c("representation", "B", "sigma2"),
    info = paste(model, "returns the documented factor schema.")
  )
  expect_identical(
    factor$representation,
    "B_sigma2_v1",
    info = paste(model, "identifies the factor representation version.")
  )
  expect_identical(
    factor$B,
    posterior_B,
    info = paste(model, "keeps B aligned with each posterior predictive draw.")
  )
  expect_identical(
    dim(factor$sigma2),
    expected_dimensions,
    info = paste(model, "returns structural variances with N by H by S dimensions.")
  )

  structural_covariance_error = 0
  for (s in seq_len(dim(posterior_B)[3])) {
    for (h in seq_len(horizon)) {
      covariance = forecasts$forecast_covariance[, , h, s]
      structural_covariance =
        factor$B[, , s] %*% covariance %*% t(factor$B[, , s])
      structural_covariance_error = max(
        structural_covariance_error,
        abs(structural_covariance - diag(factor$sigma2[, h, s]))
      )
    }
  }
  expect_true(
    structural_covariance_error < 1e-8,
    info = paste(model, "aligns B and sigma2 with every forecast covariance.")
  )
}

# for bsvar
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar$new(us_fiscal_lsuw)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2)
check_forecast_factor(ff, run_no1, 2, "PosteriorBSVAR")

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "forecast: forecast identical for normal and pipe workflow."
)

expect_true(
  is.numeric(ff$forecasts) & is.array(ff$forecasts),
  info = "forecast: returns numeric array."
)


expect_error(
  specify_bsvar$new(us_fiscal_lsuw) |> forecast(horizon = 3),
  info = "forecast: wrong input provided."
)

expect_error(
  forecast(run_no1, horizon = 1.5),
  info = "forecast: specify horizon as integer."
)




# for bsvar_msh
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar_msh$new(us_fiscal_lsuw, M = 2)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2)
check_forecast_factor(ff, run_no1, 2, "PosteriorBSVARMSH")

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar_msh$new(M = 2) |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "forecast: msh: forecast identical for normal and pipe workflow."
)

expect_identical(
  ff$forecast_mean[1,1,1], ff2$forecast_mean[1,1,1],
  info = "forecast: msh: forecast_mean identical for normal and pipe workflow."
)


expect_true(
  is.numeric(ff$forecasts) & is.array(ff$forecasts),
  info = "forecast: msh: returns numeric array."
)

expect_true(
  is.numeric(ff$forecast_mean) & is.array(ff$forecast_mean),
  info = "forecast: msh: forecast_mean: returns numeric array."
)

expect_error(
  specify_bsvar_msh$new(us_fiscal_lsuw) |> forecast(horizon = 3),
  info = "forecast: msh: wrong input provided."
)

expect_error(
  forecast(run_no1, horizon = 1.5),
  info = "forecast: msh: specify horizon as integer."
)


# for bsvar_mix
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar_mix$new(us_fiscal_lsuw, M = 2)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2)
check_forecast_factor(ff, run_no1, 2, "PosteriorBSVARMIX")

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar_mix$new(M = 2) |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "forecast: mix: forecast identical for normal and pipe workflow."
)

expect_identical(
  ff$forecast_mean[1,1,1], ff2$forecast_mean[1,1,1],
  info = "forecast: mix: forecast_mean identical for normal and pipe workflow."
)


expect_true(
  is.numeric(ff$forecasts) & is.array(ff$forecasts),
  info = "forecast: mix: returns numeric array."
)

expect_true(
  is.numeric(ff$forecast_mean) & is.array(ff$forecast_mean),
  info = "forecast: mix: forecast_mean: returns numeric array."
)

expect_error(
  specify_bsvar_msh$new(us_fiscal_lsuw) |> forecast(horizon = 3),
  info = "forecast: mix: wrong input provided."
)

expect_error(
  forecast(run_no1, horizon = 1.5),
  info = "forecast: mix: specify horizon as integer."
)


# for bsvar_sv
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar_sv$new(us_fiscal_lsuw)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2)
check_forecast_factor(ff, run_no1, 2, "PosteriorBSVARSV")

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar_sv$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "forecast: sv: forecast identical for normal and pipe workflow."
)

expect_identical(
  ff$forecast_mean[1,1,1], ff2$forecast_mean[1,1,1],
  info = "forecast: sv: forecast_mean identical for normal and pipe workflow."
)


expect_true(
  is.numeric(ff$forecasts) & is.array(ff$forecasts),
  info = "forecast: sv: returns numeric array."
)

expect_true(
  is.numeric(ff$forecast_mean) & is.array(ff$forecast_mean),
  info = "forecast: sv: forecast_mean: returns numeric array."
)

expect_error(
  specify_bsvar_msh$new(us_fiscal_lsuw) |> forecast(horizon = 3),
  info = "forecast: sv: wrong input provided."
)

expect_error(
  forecast(run_no1, horizon = 1.5),
  info = "forecast: sv: specify horizon as integer."
)


# for bsvar_sv centred
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar_sv$new(us_fiscal_lsuw, centred_sv = TRUE)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2)
check_forecast_factor(ff, run_no1, 2, "PosteriorBSVARSV (centred)")

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar_sv$new(centred_sv = TRUE) |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "forecast: sv centred: forecast identical for normal and pipe workflow."
)

expect_identical(
  ff$forecast_mean[1,1,1], ff2$forecast_mean[1,1,1],
  info = "forecast: sv centred: forecast_mean identical for normal and pipe workflow."
)






# for bsvar_t
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar_t$new(us_fiscal_lsuw)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2)
check_forecast_factor(ff, run_no1, 2, "PosteriorBSVART")

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar_t$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "forecast: t: forecast identical for normal and pipe workflow."
)

expect_identical(
  ff$forecast_mean[1,1,1], ff2$forecast_mean[1,1,1],
  info = "forecast: t: forecast_mean identical for normal and pipe workflow."
)




# conditional forecasting
################################
cf        = matrix(NA , 2, 3)
cf[,3]    = tail(us_fiscal_lsuw, 1)[3]   # conditional forecasts equal to the last gdp observation

# for bsvar
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar$new(us_fiscal_lsuw)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2, conditional_forecast = cf)
expect_false(
  "forecast_factor" %in% names(ff),
  info = "Conditional forecasts do not expose an unconditional factor."
)

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2, conditional_forecast = cf)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "conditonal forecast: forecast identical for normal and pipe workflow."
)

expect_true(
  is.numeric(ff$forecasts) & is.array(ff$forecasts),
  info = "conditonal forecast: returns numeric array."
)

expect_error(
  forecast(run_no1, horizon = 3, conditional_forecast = cf),
  info = "conditonal forecast: wrong value of horizon."
)

set.seed(20260810)
unconditional_na = forecast(
  run_no1,
  horizon = 2,
  conditional_forecast = matrix(NA_real_, 2, 3)
)
unconditional_na_seed = .Random.seed
set.seed(20260810)
unconditional_inf = forecast(
  run_no1,
  horizon = 2,
  conditional_forecast = matrix(Inf, 2, 3)
)
expect_identical(
  unconditional_inf,
  unconditional_na,
  info = "All non-finite conditioning values use the unconditional forecast path."
)
expect_identical(
  .Random.seed,
  unconditional_na_seed,
  info = "All non-finite conditioning values consume the same predictive RNG."
)

set.seed(20260810)
fully_conditional_seed = .Random.seed
fully_conditional = forecast(
  run_no1,
  horizon = 2,
  conditional_forecast = matrix(0, 2, 3)
)
expect_false(
  "forecast_factor" %in% names(fully_conditional),
  info = "Fully conditional forecasts do not expose an unconditional factor."
)
expect_identical(
  .Random.seed,
  fully_conditional_seed,
  info = "Fully conditional normal forecasts do not consume predictive-draw RNG."
)


# for bsvar_msh
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar_msh$new(us_fiscal_lsuw, M = 2)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2, conditional_forecast = cf)
expect_false(
  "forecast_factor" %in% names(ff),
  info = "Conditional forecasts do not expose an unconditional factor."
)

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar_msh$new(M = 2) |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2, conditional_forecast = cf)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "conditonal forecast: msh: forecast identical for normal and pipe workflow."
)

expect_identical(
  ff$forecast_mean[1,1,1], ff2$forecast_mean[1,1,1],
  info = "conditonal forecast: msh: forecast_mean identical for normal and pipe workflow."
)




# for bsvar_mix
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar_mix$new(us_fiscal_lsuw, M = 2)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2, conditional_forecast = cf)
expect_false(
  "forecast_factor" %in% names(ff),
  info = "Conditional forecasts do not expose an unconditional factor."
)

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar_mix$new(M = 2) |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2, conditional_forecast = cf)
)

expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "conditonal forecast: mix: forecast identical for normal and pipe workflow."
)

expect_identical(
  ff$forecast_mean[1,1,1], ff2$forecast_mean[1,1,1],
  info = "conditonal forecast: mix: forecast_mean identical for normal and pipe workflow."
)

expect_true(
  is.numeric(ff$forecasts) & is.array(ff$forecasts),
  info = "conditonal forecast: mix: returns numeric array."
)







# for bsvar_sv
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar_sv$new(us_fiscal_lsuw)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2, conditional_forecast = cf)
expect_false(
  "forecast_factor" %in% names(ff),
  info = "Conditional forecasts do not expose an unconditional factor."
)

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar_sv$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2, conditional_forecast = cf)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "conditonal forecast: sv: forecast identical for normal and pipe workflow."
)

expect_identical(
  ff$forecast_mean[1,1,1], ff2$forecast_mean[1,1,1],
  info = "conditonal forecast: sv: forecast_mean identical for normal and pipe workflow."
)




# for bsvar_t
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar_t$new(us_fiscal_lsuw)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2, conditional_forecast = cf)
expect_false(
  "forecast_factor" %in% names(ff),
  info = "Conditional forecasts do not expose an unconditional factor."
)

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar_t$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2, conditional_forecast = cf)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "conditonal forecast: t: forecast identical for normal and pipe workflow."
)

expect_identical(
  ff$forecast_mean[1,1,1], ff2$forecast_mean[1,1,1],
  info = "conditonal forecast: t: forecast_mean identical for normal and pipe workflow."
)



# for bsvar_hmsh
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar_hmsh$new(us_fiscal_lsuw)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2)
check_forecast_factor(ff, run_no1, 2, "PosteriorBSVARHMSH")
ff_conditional      <- forecast(run_no1, horizon = 2, conditional_forecast = cf)
expect_false(
  "forecast_factor" %in% names(ff_conditional),
  info = "PosteriorBSVARHMSH conditional forecasts do not expose a factor."
)

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar_hmsh$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "conditonal forecast: hmsh: forecast identical for normal and pipe workflow."
)

expect_identical(
  ff$forecast_mean[1,1,1], ff2$forecast_mean[1,1,1],
  info = "conditonal forecast: hmsh: forecast_mean identical for normal and pipe workflow."
)







# for bsvar_exh
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar_exh$new(us_fiscal_lsuw)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2)
check_forecast_factor(ff, run_no1, 2, "PosteriorBSVAREXH")
ff_conditional      <- forecast(run_no1, horizon = 2, conditional_forecast = cf)
expect_false(
  "forecast_factor" %in% names(ff_conditional),
  info = "PosteriorBSVAREXH conditional forecasts do not expose a factor."
)

set.seed(1)
suppressMessages(
  ff2              <- us_fiscal_lsuw |>
    specify_bsvar_exh$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "conditonal forecast: hmsh: forecast identical for normal and pipe workflow."
)

expect_identical(
  ff$forecast_mean[1,1,1], ff2$forecast_mean[1,1,1],
  info = "conditonal forecast: hmsh: forecast_mean identical for normal and pipe workflow."
)
