N = 1L

# These restriction rows make the prior design columns differ only through
# two 2^-128 entries.  Exact D'D contains the positive pivot 2^-255;
# sequential binary256 additions erase both half-ulp terms and report rank 1.
epsilon = 2^-128
K = 3L
T = 1L
rank_draw = .Call(
  "_bsvars_sample_A_heterosk1",
  matrix(0, N, K),
  matrix(1, N, N),
  matrix(1, 2 * N + 1L, 2L),
  matrix(1, N, T),
  matrix(0, N, T),
  matrix(0, K, T),
  list(A = matrix(0, N, K), A_V_inv = diag(K)),
  list(rbind(c(1, 0, 0), c(1, epsilon, epsilon))),
  PACKAGE = "bsvars"
)
expect_equal(dim(rank_draw), c(N, K))
expect_true(
  all(is.finite(rank_draw)),
  info = "Exact D'D retains a positive term lost by sequential binary256 sums."
)

# This is a power-of-two-scaled version of
# c(1, 2^-300, 1) %*% c(1, 1, -1) = 2^-300.  Exact D'y is 2^300 here;
# sequential binary256 accumulation produces zero after signed cancellation.
scale = 2^300
K = 1L
T = 2L
aux_A = matrix(0, N, K)
aux_B = matrix(1, N, N)
aux_hyper = matrix(1, 2 * N + 1L, 2L)
aux_sigma = matrix(1, N, T)
X = matrix(c(1, scale), K, T)
VA = list(matrix(1, 1L, K))
zero_prior = list(A = matrix(0, N, K), A_V_inv = matrix(scale^2, K, K))
signed_prior = zero_prior
signed_prior$A[] = 1

set.seed(1002)
zero_draw = .Call(
  "_bsvars_sample_A_heterosk1",
  aux_A,
  aux_B,
  aux_hyper,
  aux_sigma,
  matrix(0, N, T),
  X,
  zero_prior,
  VA,
  PACKAGE = "bsvars"
)
zero_seed = .Random.seed
set.seed(1002)
signed_draw = .Call(
  "_bsvars_sample_A_heterosk1",
  aux_A,
  aux_B,
  aux_hyper,
  aux_sigma,
  matrix(c(scale, -scale), N, T),
  X,
  signed_prior,
  VA,
  PACKAGE = "bsvars"
)
signed_seed = .Random.seed

expect_equal(
  drop((signed_draw - zero_draw) / 2^-301),
  1,
  tolerance = 1e-14,
  info = "Exact signed D'y survives strong cancellation."
)
expect_identical(
  signed_seed,
  zero_seed,
  info = "Exact dot formation consumes no random numbers."
)
