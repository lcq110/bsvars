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
  B_nu = 1L,
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
expected_norm_squared = (T + silent_prior$B_nu + 1) / 5
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
  c((orientation_nu + 1) / orientation_nu, 1 / orientation_nu)
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
