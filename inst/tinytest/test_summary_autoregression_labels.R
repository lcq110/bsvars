data(us_fiscal_lsuw)
data = us_fiscal_lsuw[1:40, 1:2]

specifications = list(
  BSVAR = function() specify_bsvar$new(data, p = 2),
  SV = function() specify_bsvar_sv$new(data, p = 2),
  EXH = function() specify_bsvar_exh$new(data, p = 2),
  MSH = function() specify_bsvar_msh$new(data, p = 2, M = 2),
  HMSH = function() specify_bsvar_hmsh$new(data, p = 2, M = 2),
  MIX = function() specify_bsvar_mix$new(data, p = 2, M = 2),
  T = function() specify_bsvar_t$new(data, p = 2)
)
expected = c("lag1_var1", "lag1_var2", "lag2_var1", "lag2_var2")

for (model in names(specifications)) {
  set.seed(361)
  specification = suppressMessages(specifications[[model]]())
  posterior = estimate(specification, S = 2, show_progress = FALSE)
  invisible(capture.output({
    posterior_summary = summary(posterior)
  }))

  expect_equal(
    rownames(posterior_summary$A[[1]])[1:4],
    expected,
    info = paste(model, "summary labels follow the lag-major coefficient order.")
  )
}
