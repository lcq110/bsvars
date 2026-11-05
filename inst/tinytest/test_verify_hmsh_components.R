N <- 1L
K <- 1L
T <- 2L
M <- 2L
draws <- 3L

Y <- matrix(c(0.5, -1), N, T)
X <- matrix(1, K, T)
xi_draw <- rbind(c(1, 0), c(0, 1))
xi_hmsh <- array(xi_draw, c(M, T, N))

posterior <- list(
  A       = array(0, c(N, K, draws)),
  B       = array(1, c(N, N, draws)),
  sigma2  = array(1, c(N, M, draws)),
  xi_cpp  = structure(rep(list(xi_hmsh), draws), dim = c(draws, 1L)),
  lambda  = array(1, c(N, T, draws))
)
prior <- list(sigma_nu = 3, sigma_s = 1)

out <- .Call(
  "_bsvars_verify_volatility_hmsh_cpp",
  posterior, prior, Y, X,
  PACKAGE = "bsvars"
)

dig2dirichlet_r <- function(x, a, b) {
  constant <- lgamma(sum(a) / 2) - sum(lgamma(a / 2)) - sum(log(b))
  kernel <- sum(0.5 * (a + 2) * (log(b) - log(x))) -
    0.5 * sum(a) * log(sum(b / x))
  constant + kernel
}

h0 <- rep(1 / M, M)
prior_nu <- rep(prior$sigma_nu, M)
prior_s <- rep(prior$sigma_s, M)
expected_denominator <- dig2dirichlet_r(h0, prior_nu, prior_s)

expect_equal(
  out$components$log_denominator,
  expected_denominator,
  tolerance = 1e-12,
  info = "HMSH prior ordinate uses the regime count in its homoskedastic point"
)

posterior_nu <- rowSums(xi_draw) + prior$sigma_nu
residuals <- as.numeric(Y)
posterior_s <- prior_s + c(residuals[1L]^2, residuals[2L]^2)
expected_numerator <- dig2dirichlet_r(h0, posterior_nu, posterior_s)

expect_equal(
  as.numeric(out$components$log_numerator_s),
  rep(expected_numerator, draws),
  tolerance = 1e-12,
  info = "HMSH posterior ordinates use the regime count in their homoskedastic point"
)
