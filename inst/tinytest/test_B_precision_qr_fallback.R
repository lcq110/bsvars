N = 2L
K = 1L
T = 1L

aux_B = diag(N)
aux_A = matrix(0, N, K)
aux_hyper = matrix(1, 2 * N + 1L, 2L)
aux_sigma = matrix(1, N, T)
Y = matrix(1e10, N, T)
X = matrix(0, K, T)
prior = list(
  B_nu = N,
  B_V_inv = diag(N)
)
VB = rep(list(diag(N)), N)

precision = prior$B_V_inv + tcrossprod(Y)
augmented_design = rbind(diag(N), t(Y))
expect_error(
  chol(precision),
  info = "The normal-equation Cholesky fails for the valid ill-conditioned input."
)
expect_true(
  min(svd(augmented_design)$d) > 0,
  info = "The augmented design remains full rank."
)

set.seed(1)
heterosk_draw = .Call(
  "_bsvars_sample_B_heterosk1",
  diag(N),
  aux_A,
  aux_hyper,
  aux_sigma,
  Y,
  X,
  prior,
  VB,
  PACKAGE = "bsvars"
)
heterosk_seed = .Random.seed
set.seed(1)
homosk_draw = .Call(
  "_bsvars_sample_B_homosk1",
  diag(N),
  aux_A,
  aux_hyper,
  Y,
  X,
  prior,
  VB,
  PACKAGE = "bsvars"
)
homosk_seed = .Random.seed

expect_equal(dim(heterosk_draw), c(N, N))
expect_equal(dim(homosk_draw), c(N, N))
expect_true(all(is.finite(heterosk_draw)))
expect_true(all(is.finite(homosk_draw)))
expect_equal(heterosk_draw, homosk_draw)
expect_identical(heterosk_seed, homosk_seed)

invalid_prior = prior
invalid_prior$B_V_inv = diag(c(1, -1))
expect_error(
  .Call(
    "_bsvars_sample_B_heterosk1",
    aux_B,
    aux_A,
    aux_hyper,
    aux_sigma,
    matrix(0, N, T),
    X,
    invalid_prior,
    VB,
    PACKAGE = "bsvars"
  )
)

silent_scale = 1e8
low_direction = c(1, -1) / sqrt(2)
high_direction = c(1, 1) / sqrt(2)
silent_Y = matrix(silent_scale, N, T)
silent_prior = prior
silent_prior$B_V_inv = 5 * diag(N)
silent_VB = list(
  matrix(low_direction, 1L, N),
  matrix(high_direction, 1L, N)
)
silent_precision = silent_prior$B_V_inv + tcrossprod(silent_Y)
projected_precision = drop(
  silent_VB[[1L]] %*% silent_precision %*% t(silent_VB[[1L]])
)
expect_true(
  projected_precision > 0 && abs(projected_precision - 5) > 1,
  info = "The projected normal equation can silently lose prior precision."
)

set.seed(4)
silent_norm_squared = replicate(5000L, {
  initial_B = rbind(low_direction, high_direction)
  sampled_B = .Call(
    "_bsvars_sample_B_heterosk1",
    initial_B,
    aux_A,
    aux_hyper,
    aux_sigma,
    silent_Y,
    X,
    silent_prior,
    silent_VB,
    PACKAGE = "bsvars"
  )
  sum(sampled_B[1L, ]^2)
})
expected_norm_squared = (
  T + silent_prior$B_nu - N + 1
) / silent_prior$B_V_inv[1L, 1L]
expect_equal(
  mean(silent_norm_squared),
  expected_norm_squared,
  tolerance = 0.04,
  info = "QR preserves the B conditional scale when inv_sympd silently succeeds."
)

orientation_N = 3L
orientation_T = 2L
orientation_restriction = diag(orientation_N)[1:2, , drop = FALSE]
orientation_VB = list(
  orientation_restriction,
  matrix(c(0, 1, 0), 1L, orientation_N),
  matrix(c(0, 0, 1), 1L, orientation_N)
)
orientation_Y = matrix(
  c(1.2, -0.3, 0.8, 0.4, 1.1, -0.9),
  orientation_N,
  orientation_T
)
orientation_sigma = matrix(
  c(0.7, 1.1, 1.4, 0.8, 1.3, 0.9),
  orientation_N,
  orientation_T
)
orientation_prior = list(
  B_nu = 3L,
  B_V_inv = matrix(
    c(4, 1, 0.5, 1, 3, 0.4, 0.5, 0.4, 2),
    orientation_N,
    orientation_N
  )
)
orientation_shocks = sweep(
  orientation_Y,
  2L,
  orientation_sigma[1L, ],
  "/"
)
orientation_precision = orientation_restriction %*%
  (orientation_prior$B_V_inv + tcrossprod(orientation_shocks)) %*%
  t(orientation_restriction)
orientation_nu = orientation_T + orientation_prior$B_nu
orientation_factor = chol(orientation_nu * solve(orientation_precision))
orientation_alpha_second_moment = diag(
  c(
    (orientation_nu - orientation_N + 1) / orientation_nu,
    1 / orientation_nu
  )
)
expected_second_moment = t(orientation_factor) %*%
  orientation_alpha_second_moment %*%
  orientation_factor

set.seed(781)
orientation_draws = replicate(10000L, {
  initial_B = diag(orientation_N)
  sampled_B = .Call(
    "_bsvars_sample_B_heterosk1",
    initial_B,
    matrix(0, orientation_N, 1L),
    matrix(1, 2 * orientation_N + 1L, 2L),
    orientation_sigma,
    orientation_Y,
    matrix(0, 1L, orientation_T),
    orientation_prior,
    orientation_VB,
    PACKAGE = "bsvars"
  )
  sampled_B[1L, 1:2]
})
empirical_second_moment = tcrossprod(orientation_draws) / ncol(orientation_draws)
expect_true(
  max(abs(empirical_second_moment - expected_second_moment)) < 0.03,
  info = paste(
    "The multivariate B draw uses the correctly oriented upper factor",
    "of the inverse precision."
  )
)

fallback_N = 2L
fallback_T = 1L
fallback_scale = 1e10
fallback_prior_precision = 5
fallback_posterior_nu = 4L
fallback_prior = list(
  B_nu = fallback_posterior_nu - fallback_T,
  B_V_inv = fallback_prior_precision * diag(fallback_N)
)
fallback_VB = list(
  diag(fallback_N),
  matrix(c(0, 1), 1L, fallback_N)
)
fallback_aux_A = matrix(0, fallback_N, 1L)
fallback_aux_hyper = matrix(1, 2 * fallback_N + 1L, 2L)
fallback_sigma = matrix(1, fallback_N, fallback_T)
fallback_X = matrix(0, 1L, fallback_T)
fallback_Y = matrix(fallback_scale, fallback_N, fallback_T)
sample_fallback_B = function(
    Y,
    draw_prior = fallback_prior,
    draw_VB = fallback_VB
) {
  .Call(
    "_bsvars_sample_B_heterosk1",
    diag(fallback_N),
    fallback_aux_A,
    fallback_aux_hyper,
    fallback_sigma,
    Y,
    fallback_X,
    draw_prior,
    draw_VB,
    PACKAGE = "bsvars"
  )
}
fallback_design = rbind(
  sqrt(fallback_prior_precision) * diag(fallback_N),
  rep(fallback_scale, fallback_N)
)
fallback_R = qr.R(qr(fallback_design, LAPACK = TRUE))
fallback_threshold = fallback_N * sqrt(.Machine$double.eps)
expect_true(
  rcond(fallback_R) < fallback_threshold / 20,
  info = "The two-dimensional fixture is well inside the fixed256 branch."
)

# Here P = v I + s^2 11'.  These closed-form expressions construct the
# canonical upper U with U' U = nu P^-1 without forming P or its inverse.
fallback_denominator = fallback_prior_precision * (
  fallback_prior_precision + 2 * fallback_scale^2
)
fallback_covariance_11 = fallback_posterior_nu * (
  fallback_prior_precision + fallback_scale^2
) / fallback_denominator
fallback_covariance_12 = -fallback_posterior_nu * fallback_scale^2 /
  fallback_denominator
fallback_covariance_determinant = fallback_posterior_nu^2 /
  fallback_denominator
fallback_U11 = sqrt(fallback_covariance_11)
fallback_U12 = fallback_covariance_12 / fallback_U11
fallback_U22 = sqrt(fallback_covariance_determinant) / fallback_U11
fallback_reference_upper = rbind(
  c(fallback_U11, fallback_U12),
  c(0, fallback_U22)
)
fallback_alpha_second_moment = diag(c(
  (fallback_posterior_nu - fallback_N + 1) / fallback_posterior_nu,
  1 / fallback_posterior_nu
))
fallback_expected_second_moment = t(fallback_reference_upper) %*%
  fallback_alpha_second_moment %*% fallback_reference_upper
fallback_transposed_second_moment = fallback_reference_upper %*%
  fallback_alpha_second_moment %*% t(fallback_reference_upper)
expect_true(
  max(abs(
    fallback_expected_second_moment - fallback_transposed_second_moment
  )) > 0.25,
  info = "The fallback fixture distinguishes upper from transposed orientation."
)

set.seed(812)
fallback_draws = replicate(6000L, {
  sampled_B = sample_fallback_B(fallback_Y)
  sampled_B[1L, ]
})
fallback_empirical_second_moment = tcrossprod(fallback_draws) /
  ncol(fallback_draws)
fallback_second_moment_error = max(abs(
  fallback_empirical_second_moment - fallback_expected_second_moment
))
expect_true(
  fallback_second_moment_error < 0.03,
  info = sprintf(
    "Fixed256 upper-factor second-moment error is %.6g.",
    fallback_second_moment_error
  )
)

safe_Y = matrix(0, fallback_N, fallback_T)
safe_design = rbind(
  sqrt(fallback_prior_precision) * diag(fallback_N),
  matrix(0, fallback_T, fallback_N)
)
safe_R = qr.R(qr(safe_design, LAPACK = TRUE))
expect_true(
  rcond(safe_R) > fallback_threshold,
  info = "The RNG comparison fixture takes the safe QR branch."
)
set.seed(913)
invisible(sample_fallback_B(safe_Y))
safe_seed = .Random.seed
set.seed(913)
invisible(sample_fallback_B(fallback_Y))
fallback_seed = .Random.seed
expect_identical(
  fallback_seed,
  safe_seed,
  info = "The fixed256 factorization consumes no random numbers."
)

dependent_prior = list(
  B_nu = 3L,
  B_V_inv = matrix(c(2, 0.25, 0.25, 1), fallback_N, fallback_N)
)
expect_silent(
  chol(dependent_prior$B_V_inv),
  info = "The dependent-restriction fixture has a full-SPD prior."
)
dependent_restriction = rbind(c(1, 0), c(2, 0))
dependent_VB = list(
  dependent_restriction,
  matrix(c(0, 1), 1L, fallback_N)
)
expect_true(
  qr(dependent_restriction)$rank < nrow(dependent_restriction),
  info = "The restriction rows are linearly dependent."
)
expect_error(
  sample_fallback_B(safe_Y, dependent_prior, dependent_VB),
  pattern = "multiprecision Cholesky is not positive definite",
  info = "A singular restricted precision errors instead of being repaired."
)
