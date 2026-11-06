data(us_fiscal_lsuw)
data = us_fiscal_lsuw[1:40, 1:2]

specifications = list(
  BSVAR = function() specify_bsvar$new(data),
  T = function() specify_bsvar_t$new(data),
  SV = function() specify_bsvar_sv$new(data),
  MSH = function() specify_bsvar_msh$new(data, M = 2),
  HMSH = function() specify_bsvar_hmsh$new(data, M = 2),
  EXH = function() specify_bsvar_exh$new(
    data,
    variance_regimes = rep(1:2, length.out = nrow(data))
  )
)

for (model in names(specifications)) {
  specification = suppressMessages(specifications[[model]]())
  K = ncol(specification$prior$A_V_inv)
  specification$prior$A_V_inv[K, K] = -1e12

  expect_error(
    estimate(specification, S = 2, show_progress = FALSE),
    pattern = "chol",
    info = paste(model, "propagates a failed conditional parameter draw.")
  )
}
