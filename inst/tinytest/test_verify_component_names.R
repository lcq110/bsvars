N = 1L
T = 3L
S = 2L
Y = matrix(c(0.5, -1, 1.5), N, T)
X = matrix(1, 1, T)
prior = list(sv_a_ = 1, sv_s_ = 1)

posterior = list(
  posterior = list(
    B = array(1, c(N, N, S)),
    A = array(0, c(N, 1, S)),
    h = array(0, c(N, T, S)),
    S = array(0, c(N, T, S)),
    sigma2_omega = matrix(1, N, S),
    s_ = matrix(1, N, S),
    lambda = array(1, c(N, T, S))
  ),
  last_draw = list(
    centred_sv = FALSE,
    prior = list(get_prior = function() prior),
    data_matrices = list(Y = Y, X = X)
  )
)
class(posterior) = "PosteriorBSVARSV"

volatility = verify_volatility(posterior)
identification = verify_identification(posterior)

expect_identical(
  names(volatility)[1:2],
  c("logSDDR", "logSDDR_se"),
  info = "Volatility verification exposes the documented numerical-error component name."
)
expect_identical(
  names(identification),
  c("logSDDR", "logSDDR_se"),
  info = "Identification verification exposes the documented numerical-error component name."
)
