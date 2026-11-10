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
  hypothesis_second, posterior_homosk, prior_scaled, VA_all_free, Y0, X0,
  PACKAGE = "bsvars"
)
set.seed(21)
out_identity <- .Call(
  "_bsvars_verify_autoregressive_homosk_cpp",
  hypothesis_second, posterior_homosk, prior_identity, VA_all_free, Y0, X0,
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

# The public method must forward VA to the native routine.
posterior_object <- list(
  posterior = posterior_homosk,
  last_draw = list(
    starting_values = list(A = matrix(0, N, K)),
    identification  = list(VA = VA_all_free),
    prior           = list(get_prior = function() prior_scaled),
    data_matrices   = list(Y = Y0, X = X0)
  )
)
class(posterior_object) <- "PosteriorBSVAR"

set.seed(22)
out_public <- verify_autoregression(
  posterior_object,
  matrix(c(NA, 0), N, K)
)
set.seed(22)
out_native <- .Call(
  "_bsvars_verify_autoregressive_homosk_cpp",
  hypothesis_second, posterior_homosk, prior_scaled, VA_all_free, Y0, X0,
  PACKAGE = "bsvars"
)

expect_equal(
  out_public$components$log_denominator_s,
  out_native$components$log_denominator_s,
  tolerance = 1e-12,
  info = "public autoregression verification forwards free-coordinate restrictions"
)

# A fixed regressor cannot alter the ordinate of a tested free coefficient.
VA_restricted <- list(matrix(c(1, 0), 1L, K))
hypothesis_first <- matrix(c(0, 999), N, K)
x1 <- c(1, 2, -1, 0.5)
X_without_fixed <- rbind(x1, rep(0, T))
X_with_fixed    <- rbind(x1, x1)
Y <- matrix(c(0.5, -1, 1.5, -0.75), N, T)

posterior_heterosk <- posterior_homosk
posterior_heterosk$sigma  <- array(1, c(N, T, draws))
posterior_heterosk$lambda <- array(1, c(N, T, draws))

prior_restricted <- prior_identity
prior_restricted$A_V_inv <- matrix(c(2, 1, 1, 2), K, K)

out_homosk_without_fixed <- .Call(
  "_bsvars_verify_autoregressive_homosk_cpp",
  hypothesis_first, posterior_homosk, prior_restricted, VA_restricted,
  Y, X_without_fixed,
  PACKAGE = "bsvars"
)
out_homosk_with_fixed <- .Call(
  "_bsvars_verify_autoregressive_homosk_cpp",
  hypothesis_first, posterior_homosk, prior_restricted, VA_restricted,
  Y, X_with_fixed,
  PACKAGE = "bsvars"
)

expect_equal(
  out_homosk_without_fixed$components$log_numerator_s,
  out_homosk_with_fixed$components$log_numerator_s,
  tolerance = 1e-12,
  info = "homoskedastic posterior ordinate is computed in free coordinates"
)

out_without_fixed <- .Call(
  "_bsvars_verify_autoregressive_heterosk_cpp",
  hypothesis_first, posterior_heterosk, prior_restricted, VA_restricted,
  Y, X_without_fixed,
  PACKAGE = "bsvars"
)
out_with_fixed <- .Call(
  "_bsvars_verify_autoregressive_heterosk_cpp",
  hypothesis_first, posterior_heterosk, prior_restricted, VA_restricted,
  Y, X_with_fixed,
  PACKAGE = "bsvars"
)

expect_equal(
  out_without_fixed$components$log_numerator_s,
  out_with_fixed$components$log_numerator_s,
  tolerance = 1e-12,
  info = "heteroskedastic posterior ordinate is computed in free coordinates"
)
