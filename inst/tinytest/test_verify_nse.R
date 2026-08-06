N <- 1L
K <- 1L
T <- 2L
M <- 2L

Y <- matrix(c(0.5, -1), N, T)
X <- matrix(1, K, T)

make_common_posterior <- function(draws) {
  list(
    A = array(0, c(N, K, draws)),
    B = array(1, c(N, N, draws))
  )
}

make_sv_posterior <- function(draws) {
  c(
    make_common_posterior(draws),
    list(
      h            = array(rep(c(0.2, -0.1), draws), c(N, T, draws)),
      S            = array(0, c(N, T, draws)),
      sigma2_omega = matrix(1, N, draws),
      s_           = matrix(1, N, draws),
      lambda       = array(1, c(N, T, draws))
    )
  )
}

xi_draw <- rbind(c(1, 0), c(0, 1))

make_msh_posterior <- function(draws) {
  c(
    make_common_posterior(draws),
    list(
      sigma2 = array(1, c(N, M, draws)),
      xi     = array(rep(xi_draw, draws), c(M, T, draws)),
      lambda = array(1, c(N, T, draws))
    )
  )
}

make_hmsh_posterior <- function(draws) {
  xi_hmsh <- array(xi_draw, c(M, T, N))
  c(
    make_common_posterior(draws),
    list(
      sigma2 = array(1, c(N, M, draws)),
      xi_cpp = structure(rep(list(xi_hmsh), draws), dim = c(draws, 1L)),
      lambda = array(1, c(N, T, draws))
    )
  )
}

make_ar_posterior <- function(draws) {
  c(
    make_common_posterior(draws),
    list(
      hyper  = array(1, c(2 * N + 1L, 2L, draws)),
      sigma  = array(1, c(N, T, draws)),
      lambda = array(1, c(N, T, draws))
    )
  )
}

prior_sv <- list(sv_a_ = 1, sv_s_ = 1)
prior_msh <- list(sigma_nu = 3, sigma_s = 1)
prior_ar <- list(
  hyper_nu_A  = 5,
  hyper_a_A   = 1,
  hyper_s_AA  = 1,
  hyper_nu_AA = 5,
  A           = matrix(0, N, K),
  A_V_inv     = matrix(1, K, K)
)
hypothesis <- matrix(0, N, K)
VA <- list(matrix(1, N, K))

short_draws <- 2L

out_sv_short <- .Call(
  "_bsvars_verify_volatility_sv_cpp",
  make_sv_posterior(short_draws), prior_sv, Y, X, FALSE,
  PACKAGE = "bsvars"
)
out_msh_short <- .Call(
  "_bsvars_verify_volatility_msh_cpp",
  make_msh_posterior(short_draws), prior_msh, Y, X,
  PACKAGE = "bsvars"
)
out_hmsh_short <- .Call(
  "_bsvars_verify_volatility_hmsh_cpp",
  make_hmsh_posterior(short_draws), prior_msh, Y, X,
  PACKAGE = "bsvars"
)
set.seed(301)
out_ar_homosk_short <- .Call(
  "_bsvars_verify_autoregressive_homosk_cpp",
  hypothesis, make_ar_posterior(short_draws), prior_ar, VA, Y, X,
  PACKAGE = "bsvars"
)
set.seed(302)
out_ar_heterosk_short <- .Call(
  "_bsvars_verify_autoregressive_heterosk_cpp",
  hypothesis, make_ar_posterior(short_draws), prior_ar, VA, Y, X,
  PACKAGE = "bsvars"
)

short_outputs <- list(
  sv         = list(out_sv_short$logSDDR_se, out_sv_short$components$se_components),
  msh        = list(out_msh_short$logSDDR_se, out_msh_short$components$se_components),
  hmsh       = list(out_hmsh_short$logSDDR_se, out_hmsh_short$components$se_components),
  ar_homosk  = list(out_ar_homosk_short$log_SDDR_se, out_ar_homosk_short$components$se_components),
  ar_heterosk = list(out_ar_heterosk_short$log_SDDR_se, out_ar_heterosk_short$components$se_components)
)

for (path in names(short_outputs)) {
  expect_true(
    all(is.nan(short_outputs[[path]][[1L]])),
    info = sprintf("%s NSE is NaN when fewer than 60 draws are available", path)
  )
  expect_true(
    all(is.nan(short_outputs[[path]][[2L]])),
    info = sprintf("%s NSE components are NaN when fewer than 60 draws are available", path)
  )
}

log_mean_r <- function(x) {
  anchor <- max(x)
  anchor + log(mean(exp(x - anchor)))
}

batch_components <- function(out) {
  vapply(seq_len(30L), function(batch) {
    indices <- (2L * batch - 1L):(2L * batch)
    log_mean_r(out$components$log_numerator_s[1L, indices]) -
      log_mean_r(out$components$log_denominator_s[1L, indices])
  }, numeric(1L))
}

long_draws <- 60L

set.seed(303)
out_ar_homosk <- .Call(
  "_bsvars_verify_autoregressive_homosk_cpp",
  hypothesis, make_ar_posterior(long_draws), prior_ar, VA, Y, X,
  PACKAGE = "bsvars"
)
expected_ar_homosk <- batch_components(out_ar_homosk)

expect_equal(
  as.numeric(out_ar_homosk$components$se_components),
  expected_ar_homosk,
  tolerance = 1e-12,
  info = "homoskedastic autoregression NSE components start from zero"
)
expect_equal(
  out_ar_homosk$log_SDDR_se,
  sqrt(mean((expected_ar_homosk - mean(expected_ar_homosk))^2)),
  tolerance = 1e-12,
  info = "homoskedastic autoregression NSE uses the initialized components"
)

set.seed(304)
out_ar_heterosk <- .Call(
  "_bsvars_verify_autoregressive_heterosk_cpp",
  hypothesis, make_ar_posterior(long_draws), prior_ar, VA, Y, X,
  PACKAGE = "bsvars"
)
expected_ar_heterosk <- batch_components(out_ar_heterosk)

expect_equal(
  as.numeric(out_ar_heterosk$components$se_components),
  expected_ar_heterosk,
  tolerance = 1e-12,
  info = "heteroskedastic autoregression NSE components start from zero"
)
expect_equal(
  out_ar_heterosk$log_SDDR_se,
  sqrt(mean((expected_ar_heterosk - mean(expected_ar_heterosk))^2)),
  tolerance = 1e-12,
  info = "heteroskedastic autoregression NSE uses the initialized components"
)

set.seed(305)
out_sv <- .Call(
  "_bsvars_verify_volatility_sv_cpp",
  make_sv_posterior(long_draws), prior_sv, Y, X, TRUE,
  PACKAGE = "bsvars"
)
prior_sv_scaled <- prior_sv
prior_sv_scaled$sv_s_ <- 4
set.seed(305)
out_sv_scaled <- .Call(
  "_bsvars_verify_volatility_sv_cpp",
  make_sv_posterior(long_draws), prior_sv_scaled, Y, X, TRUE,
  PACKAGE = "bsvars"
)

expect_equal(
  as.numeric(
    out_sv_scaled$components$se_components -
      out_sv$components$se_components
  ),
  rep(log(2), 30L),
  tolerance = 1e-12,
  info = "SV NSE batches subset the original prior-scale draws"
)
