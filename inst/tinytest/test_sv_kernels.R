.sv_mixture_cdf <- function(x) {
  .Call(bsvars:::`_bsvars_find_mixture_indicator_cdf`, x)
}

.sv_mixture_indicators <- function(cdf, T) {
  .Call(
    bsvars:::`_bsvars_inverse_transform_sampling`,
    cdf,
    as.integer(T)
  )
}

.centred_sv_draw <- function(
    seed,
    h,
    u,
    prior = list(sv_a_ = 1, sv_s_ = 0.1),
    s_ = 0.05,
    sample_s_ = TRUE
) {
  set.seed(seed)
  .Call(
    bsvars:::`_bsvars_svar_ce1`,
    h + 0,
    0.5,
    2,
    0.1,
    1,
    s_,
    rep.int(0L, length(h)),
    u,
    prior,
    sample_s_
  )
}


# B04: mixture CDF uses the normalized Gaussian component densities.
alpha <- c(
  1.92677, 1.34744, 0.73504, 0.02266, -0.85173,
  -1.97278, -3.46788, -5.55246, -8.68384, -14.65000
)
variance <- c(
  0.11265, 0.17788, 0.26768, 0.40611, 0.62699,
  0.98583, 1.57469, 2.54498, 4.16591, 7.33342
)
probability <- c(
  0.00609, 0.04775, 0.13057, 0.20674, 0.22715,
  0.18842, 0.12047, 0.05591, 0.01575, 0.00115
)
data_normal <- c(-3, 0, 2)

expected_cdf <- vapply(
  data_normal,
  function(x) {
    log_weight <- log(probability) - 0.5 * log(variance) -
      0.5 * (x - alpha)^2 / variance
    weight <- exp(log_weight - max(log_weight))
    cumsum(weight / sum(weight))
  },
  numeric(10)
)

expect_equal(
  as.numeric(.sv_mixture_cdf(data_normal)),
  as.numeric(expected_cdf),
  tolerance = 1e-12,
  info = "find_mixture_indicator_cdf: CDF matches normalized Gaussian mixture weights."
)
