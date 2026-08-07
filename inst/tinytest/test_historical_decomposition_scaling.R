N = 2
T = 2
p = 1
B = diag(c(2, 3))
A = matrix(0, N, N * p + 1)
Y = matrix(1, N, T)
X = rbind(matrix(0, N, T), 1)

posterior = list(
  posterior = list(
    B = array(B, c(N, N, 1)),
    A = array(A, c(N, ncol(A), 1))
  ),
  last_draw = list(
    data_matrices = list(Y = Y, X = X),
    p = p
  )
)

classes = c(
  "PosteriorBSVAR",
  "PosteriorBSVAREXH",
  "PosteriorBSVARMSH",
  "PosteriorBSVARHMSH",
  "PosteriorBSVARMIX",
  "PosteriorBSVARSV",
  "PosteriorBSVART"
)

for (model_class in classes) {
  class(posterior) = model_class
  decomposition = compute_historical_decompositions(
    posterior,
    show_progress = FALSE
  )
  reconstructed_residuals = vapply(
    1:T,
    function(t) rowSums(decomposition[, , t, 1]),
    numeric(N)
  )

  expect_equal(
    reconstructed_residuals,
    Y,
    info = paste(model_class, "uses unstandardized impulse responses in historical decompositions.")
  )
}
