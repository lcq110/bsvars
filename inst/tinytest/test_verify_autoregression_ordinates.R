N     <- 1L
K     <- 2L
T     <- 4L
draws <- 2L

B <- array(1, c(N, N, draws))
A <- array(0, c(N, K, draws))
posterior_homosk <- list(
  A     = A,
  B     = B,
  hyper = array(1, c(2 * N + 1L, 2L, draws))
)

prior_identity <- list(
  hyper_nu_A  = 5,
  hyper_a_A   = 1,
  hyper_s_AA  = 1,
  hyper_nu_AA = 5,
  A           = matrix(0, N, K),
  A_V_inv     = diag(K)
)

# The prior ordinate must use A_V_inv rather than an identity covariance.
prior_scaled <- prior_identity
prior_scaled$A_V_inv <- diag(c(1, 4))
hypothesis_second <- matrix(c(999, 0), N, K)
VA_all_free <- list(diag(K))
Y0 <- matrix(0, N, T)
X0 <- matrix(0, K, T)

set.seed(21)
out_scaled <- .Call(
  "_bsvars_verify_autoregressive_homosk_cpp",
  hypothesis_second, posterior_homosk, prior_scaled, Y0, X0,
  PACKAGE = "bsvars"
)
set.seed(21)
out_identity <- .Call(
  "_bsvars_verify_autoregressive_homosk_cpp",
  hypothesis_second, posterior_homosk, prior_identity, Y0, X0,
  PACKAGE = "bsvars"
)

expect_equal(
  as.numeric(
    out_scaled$components$log_denominator_s -
      out_identity$components$log_denominator_s
  ),
  rep(log(2), draws),
  tolerance = 1e-12,
  info = "autoregression prior ordinate uses the configured precision"
)
