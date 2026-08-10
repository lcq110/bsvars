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
  info = "The QR fallback draws from a valid ill-conditioned conditional posterior."
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
  info = "The fallback uses the correct conditional-posterior location."
)
