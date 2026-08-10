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

expect_equal(dim(heterosk_draw), c(N, N))
expect_equal(dim(homosk_draw), c(N, N))
expect_true(all(is.finite(heterosk_draw)))
expect_true(all(is.finite(homosk_draw)))
expect_equal(heterosk_draw, homosk_draw)

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
