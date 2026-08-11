# With no observations and the default B_nu = N, the generalized-normal
# prior reduces to a normal prior.  In particular, the first diagonal element
# of a lower-triangular B has unit second moment.
prior_only_N = 4L
prior_only_B = diag(prior_only_N)
prior_only_A = matrix(0, prior_only_N, 1L)
prior_only_hyper = matrix(1, 2 * prior_only_N + 1L, 2L)
prior_only_Y = matrix(numeric(), prior_only_N, 0L)
prior_only_X = matrix(numeric(), 1L, 0L)
prior_only_prior = list(
  B_nu = prior_only_N,
  B_V_inv = diag(prior_only_N)
)
prior_only_VB = lapply(seq_len(prior_only_N), function(n) {
  diag(prior_only_N)[seq_len(n), , drop = FALSE]
})

set.seed(1701)
prior_only_b11_squared = replicate(4000L, {
  sampled_B = .Call(
    "_bsvars_sample_B_homosk1",
    prior_only_B,
    prior_only_A,
    prior_only_hyper,
    prior_only_Y,
    prior_only_X,
    prior_only_prior,
    prior_only_VB,
    PACKAGE = "bsvars"
  )
  sampled_B[1L, 1L]^2
})
expect_equal(
  mean(prior_only_b11_squared),
  1,
  tolerance = 0.08,
  info = paste(
    "The default generalized-normal B prior retains its documented",
    "normal-prior boundary."
  )
)

# For B_nu > N, the same determinant power also enters the conditional
# degrees of freedom of the row-specific B shrinkage parameter.
hyper_N = 2L
hyper_K = 2L
hyper_B = matrix(c(1.2, 0, 0.3, 0.9), hyper_N, hyper_N, byrow = TRUE)
hyper_A = matrix(c(0.5, -0.2, 0.1, 0.7), hyper_N, hyper_K, byrow = TRUE)
hyper_VB = list(
  matrix(c(1, 0), 1L, hyper_N),
  diag(hyper_N)
)
hyper_VA = rep(list(diag(hyper_K)), hyper_N)
hyper_template = matrix(
  c(
    1.1, 1.3, 0.8, 0.9, 1.5,
    1.2, 1.4, 0.7, 1.0, 1.6
  ),
  2 * hyper_N + 1L,
  2L
)
hyper_prior = list(
  B_nu = hyper_N + 3L,
  B_V_inv = matrix(c(2, 0.25, 0.25, 1.5), hyper_N, hyper_N),
  A = matrix(c(0.1, 0, 0, -0.1), hyper_N, hyper_K, byrow = TRUE),
  A_V_inv = matrix(c(1.5, 0.2, 0.2, 1.2), hyper_K, hyper_K),
  hyper_nu_B = 7,
  hyper_a_B = 2.5,
  hyper_s_BB = 3,
  hyper_nu_BB = 4,
  hyper_nu_A = 6,
  hyper_a_A = 2,
  hyper_s_AA = 2.5,
  hyper_nu_AA = 5
)

replay_hyperparameters = function(aux_hyper) {
  output = aux_hyper

  scale = hyper_prior$hyper_s_BB + 2 * sum(output[(hyper_N + 1L):(2 * hyper_N), 1L])
  shape = hyper_prior$hyper_nu_BB + 2 * hyper_N * hyper_prior$hyper_a_B
  output[2 * hyper_N + 1L, 1L] = scale / rchisq(1L, shape)

  scale = hyper_prior$hyper_s_AA + 2 * sum(output[(hyper_N + 1L):(2 * hyper_N), 2L])
  shape = hyper_prior$hyper_nu_AA + 2 * hyper_N * hyper_prior$hyper_a_A
  output[2 * hyper_N + 1L, 2L] = scale / rchisq(1L, shape)

  for (n in seq_len(hyper_N)) {
    scale = 1 / (
      1 / (2 * output[n, 1L]) + 1 / output[2 * hyper_N + 1L, 1L]
    )
    shape = hyper_prior$hyper_a_B + hyper_prior$hyper_nu_B / 2
    output[hyper_N + n, 1L] = rgamma(1L, shape, scale = scale)

    quadratic = drop(
      hyper_B[n, , drop = FALSE] %*%
        hyper_prior$B_V_inv %*%
        t(hyper_B[n, , drop = FALSE])
    )
    scale = output[hyper_N + n, 1L] + quadratic
    shape = hyper_prior$hyper_nu_B + nrow(hyper_VB[[n]]) +
      hyper_prior$B_nu - hyper_N
    output[n, 1L] = scale / rchisq(1L, shape)

    scale = 1 / (
      1 / (2 * output[n, 2L]) + 1 / output[2 * hyper_N + 1L, 2L]
    )
    shape = hyper_prior$hyper_a_A + hyper_prior$hyper_nu_A / 2
    output[hyper_N + n, 2L] = rgamma(1L, shape, scale = scale)

    difference = hyper_A[n, , drop = FALSE] -
      hyper_prior$A[n, , drop = FALSE]
    quadratic = drop(
      difference %*% hyper_prior$A_V_inv %*% t(difference)
    )
    scale = output[hyper_N + n, 2L] + quadratic
    shape = hyper_prior$hyper_nu_A + nrow(hyper_VA[[n]])
    output[n, 2L] = scale / rchisq(1L, shape)
  }

  output
}

set.seed(2701)
expected_hyper = replay_hyperparameters(hyper_template)
set.seed(2701)
sampled_hyper = .Call(
  "_bsvars_sample_hyperparameters",
  hyper_template + 0,
  hyper_B,
  hyper_A,
  hyper_VB,
  hyper_VA,
  hyper_prior,
  PACKAGE = "bsvars"
)
expect_equal(
  sampled_hyper,
  expected_hyper,
  tolerance = 1e-14,
  info = paste(
    "The B shrinkage conditional uses the same determinant power as the",
    "generalized-normal B conditional."
  )
)
