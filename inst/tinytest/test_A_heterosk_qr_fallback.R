N = 1L
K = 2L
T = 2L

aux_A = matrix(0, N, K)
aux_B = matrix(1, N, N)
aux_hyper = matrix(1, 2 * N + 1, 2)
aux_sigma = matrix(c(1e-10, 1), N, T)
Y = matrix(0, N, T)
X = matrix(c(1, 1, 0, 1), K, T)
prior = list(
  A = matrix(0, N, K),
  A_V_inv = diag(K)
)
VA = list(diag(K))

weighted_design = t(X) / drop(aux_sigma)
precision = prior$A_V_inv + crossprod(weighted_design)
augmented_design = rbind(diag(K), weighted_design)
expect_error(
  chol(precision),
  info = "The normal-equation Cholesky fails for the valid ill-conditioned input."
)
expect_true(
  min(svd(augmented_design)$d) > 0,
  info = "The augmented design remains full rank."
)

set.seed(1)
draw = .Call(
  "_bsvars_sample_A_heterosk1",
  aux_A,
  aux_B,
  aux_hyper,
  aux_sigma,
  Y,
  X,
  prior,
  VA,
  PACKAGE = "bsvars"
)

expect_equal(dim(draw), c(N, K))
expect_true(
  all(is.finite(draw)),
  info = "QR draws from a valid ill-conditioned conditional posterior."
)

shifted_Y = matrix(c(1, -2), N, T)
shifted_prior = prior
shifted_prior$A = matrix(c(0.5, -0.25), N, K)
augmented_response = c(
  drop(shifted_prior$A),
  drop(shifted_Y) / drop(aux_sigma)
)
expected_mean_shift = qr.coef(
  qr(augmented_design, LAPACK = TRUE),
  augmented_response
)

set.seed(2)
zero_draw = .Call(
  "_bsvars_sample_A_heterosk1",
  aux_A,
  aux_B,
  aux_hyper,
  aux_sigma,
  Y,
  X,
  prior,
  VA,
  PACKAGE = "bsvars"
)
set.seed(2)
shifted_draw = .Call(
  "_bsvars_sample_A_heterosk1",
  aux_A,
  aux_B,
  aux_hyper,
  aux_sigma,
  shifted_Y,
  X,
  shifted_prior,
  VA,
  PACKAGE = "bsvars"
)

expect_equal(
  drop(shifted_draw - zero_draw),
  drop(expected_mean_shift),
  tolerance = 1e-8,
  info = "QR uses the correct conditional-posterior location."
)

silent_scale = 1e8
silent_prior = prior
silent_prior$A_V_inv = 5 * diag(K)
silent_X = matrix(c(silent_scale, silent_scale), K, 1L)
silent_sigma = matrix(1, N, 1L)
silent_Y = matrix(0, N, 1L)
silent_precision = silent_prior$A_V_inv + tcrossprod(drop(silent_X))
expect_silent(
  chol(silent_precision),
  info = "The inaccurate normal-equation Cholesky can report success."
)

target_mean = c(1, -1) / sqrt(2)
shifted_silent_prior = silent_prior
shifted_silent_prior$A = matrix(target_mean, N, K)
set.seed(3)
zero_silent_draw = .Call(
  "_bsvars_sample_A_heterosk1",
  aux_A,
  aux_B,
  aux_hyper,
  silent_sigma,
  silent_Y,
  silent_X,
  silent_prior,
  VA,
  PACKAGE = "bsvars"
)
set.seed(3)
shifted_silent_draw = .Call(
  "_bsvars_sample_A_heterosk1",
  aux_A,
  aux_B,
  aux_hyper,
  silent_sigma,
  silent_Y,
  silent_X,
  shifted_silent_prior,
  VA,
  PACKAGE = "bsvars"
)
expect_equal(
  drop(shifted_silent_draw - zero_silent_draw),
  target_mean,
  tolerance = 1e-7,
  info = "QR remains accurate when the normal-equation Cholesky silently succeeds."
)

set.seed(5)
heterosk_silent_draw = .Call(
  "_bsvars_sample_A_heterosk1",
  aux_A,
  aux_B,
  aux_hyper,
  silent_sigma,
  silent_Y,
  silent_X,
  shifted_silent_prior,
  VA,
  PACKAGE = "bsvars"
)
heterosk_seed = .Random.seed
set.seed(5)
homosk_silent_draw = .Call(
  "_bsvars_sample_A_homosk1",
  aux_A,
  aux_B,
  aux_hyper,
  silent_Y,
  silent_X,
  shifted_silent_prior,
  VA,
  PACKAGE = "bsvars"
)
homosk_seed = .Random.seed
expect_equal(heterosk_silent_draw, homosk_silent_draw)
expect_identical(heterosk_seed, homosk_seed)

covariance_K = 3L
covariance_pattern = c(TRUE, TRUE, FALSE)
covariance_aux_A = matrix(0, N, covariance_K)
covariance_aux_hyper = matrix(1, 2 * N + 1L, 2L)
covariance_sigma = matrix(1, N, 1L)
covariance_Y = matrix(0, N, 1L)
covariance_X = matrix(0, covariance_K, 1L)
covariance_prior = list(
  A = matrix(0, N, covariance_K),
  A_V_inv = matrix(
    c(
      4, 1.5, 0,
      1.5, 1, 0,
      0, 0, 2
    ),
    covariance_K,
    covariance_K
  )
)
covariance_VA = list(matrix(
  diag(covariance_K)[covariance_pattern, , drop = FALSE],
  ncol = covariance_K
))
restricted_precision = covariance_VA[[1L]] %*%
  covariance_prior$A_V_inv %*% t(covariance_VA[[1L]])
theoretical_covariance = solve(restricted_precision)
precision_chol = chol(restricted_precision)
wrong_orientation_covariance = solve(t(precision_chol)) %*%
  solve(precision_chol)
orientation_gap = max(abs(
  wrong_orientation_covariance - theoretical_covariance
))
expect_true(
  orientation_gap > 0.5,
  info = "The fixture distinguishes R^-1 from R^-T noise orientation."
)

set.seed(6)
covariance_draws = replicate(50000L, {
  sampled_A = .Call(
    "_bsvars_sample_A_heterosk1",
    covariance_aux_A,
    aux_B,
    covariance_aux_hyper,
    covariance_sigma,
    covariance_Y,
    covariance_X,
    covariance_prior,
    covariance_VA,
    PACKAGE = "bsvars"
  )
  drop(sampled_A[1L, covariance_pattern])
})
empirical_covariance = cov(t(covariance_draws))
covariance_error = max(abs(empirical_covariance - theoretical_covariance))
expect_true(
  covariance_error < 0.06,
  info = sprintf(
    "Empirical covariance differs from the QR target by %.6g.",
    covariance_error
  )
)
