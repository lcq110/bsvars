# HMSH FEVDs must use the terminal state from the matching posterior draw.
N = 2L
M = 2L
S = 2L
T = 1L
p = 1L

impact = matrix(c(1, 0, 1, 1), N, N)
B      = solve(impact)

posterior_B = array(rep(B, S), c(N, N, S))
posterior_A = array(0, c(N, N * p + 1, S))
for (s in 1:S) posterior_A[, 1:N, s] = 0.5 * diag(N)

transition = array(0, c(M, M, N))
for (n in 1:N) transition[, , n] = diag(M)
posterior_PR_TR = matrix(
  replicate(S, transition, simplify = FALSE),
  S,
  1
)

xi_draw_1 = array(0, c(M, T, N))
xi_draw_1[, 1, 1] = c(1, 0)
xi_draw_1[, 1, 2] = c(0, 1)
xi_draw_2 = array(0, c(M, T, N))
xi_draw_2[, 1, 1] = c(0, 1)
xi_draw_2[, 1, 2] = c(1, 0)

posterior_xi_cpp = matrix(list(xi_draw_1, xi_draw_2), S, 1)
posterior_xi     = array(0, c(M, T, N, S))
posterior_xi[, , , 1] = xi_draw_1
posterior_xi[, , , 2] = xi_draw_2

posterior_sigma2 = array(NA_real_, c(N, M, S))
posterior_sigma2[, , 1] = matrix(c(1, 4, 1, 4), N, M, byrow = TRUE)
posterior_sigma2[, , 2] = matrix(c(1, 4, 1, 4), N, M, byrow = TRUE)
posterior_sigma = array(NA_real_, c(N, T, S))
posterior_sigma[, 1, 1] = c(1, 2)
posterior_sigma[, 1, 2] = c(2, 1)

hmsh_posterior = list(
  posterior = list(
    B         = posterior_B,
    A         = posterior_A,
    PR_TR_cpp = posterior_PR_TR,
    sigma2    = posterior_sigma2,
    xi_cpp    = posterior_xi_cpp,
    xi        = posterior_xi,
    sigma     = posterior_sigma,
    lambda    = array(1, c(N, T, S)),
    df        = matrix(10, N, S)
  ),
  last_draw = list(
    p             = p,
    data_matrices = list(
      Y = matrix(0, N, T, dimnames = list(c("y1", "y2"), "t1"))
    ),
    get_normal = function() TRUE
  )
)
class(hmsh_posterior) = "PosteriorBSVARHMSH"

hmsh_fevd = compute_variance_decompositions(hmsh_posterior, horizon = 1)

expect_equal(
  unname(hmsh_fevd[1, , 2, 1]),
  c(20, 80),
  tolerance = 1e-10,
  info = "HMSH FEVD: draw 1 uses draw 1's terminal states."
)
expect_equal(
  unname(hmsh_fevd[1, , 2, 2]),
  c(80, 20),
  tolerance = 1e-10,
  info = "HMSH FEVD: draw 2 uses draw 2's terminal states."
)


# Student-t HMSH impact FEVDs must include the terminal latent scale.
student_hmsh_posterior = hmsh_posterior
student_hmsh_posterior$posterior$sigma2[] = 1
student_hmsh_posterior$posterior$sigma[]  = 1
student_hmsh_posterior$posterior$lambda[, 1, 1] = c(4, 1)
student_hmsh_posterior$posterior$lambda[, 1, 2] = c(1, 4)
student_hmsh_posterior$last_draw$get_normal = function() FALSE

set.seed(20260806)
student_hmsh_fevd = compute_variance_decompositions(
  student_hmsh_posterior,
  horizon = 1
)

expect_equal(
  unname(student_hmsh_fevd[1, , 1, 1]),
  c(80, 20),
  tolerance = 1e-10,
  info = "Student-t HMSH FEVD: draw 1 includes terminal lambda."
)
expect_equal(
  unname(student_hmsh_fevd[1, , 1, 2]),
  c(20, 80),
  tolerance = 1e-10,
  info = "Student-t HMSH FEVD: draw 2 includes terminal lambda."
)
