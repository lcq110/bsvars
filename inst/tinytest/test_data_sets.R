data(us_fiscal_sww)

expect_true(
  inherits(us_fiscal_sww, "ts"),
  info = "us_fiscal_sww is a time series object."
)
expect_identical(
  dim(us_fiscal_sww),
  c(266L, 10L),
  info = "us_fiscal_sww has the expected dimensions."
)
expect_equal(
  frequency(us_fiscal_sww),
  4,
  info = "us_fiscal_sww is quarterly."
)
expect_identical(
  as.integer(start(us_fiscal_sww)),
  c(1959L, 2L),
  info = "us_fiscal_sww starts in 1959 Q2."
)
expect_identical(
  as.integer(end(us_fiscal_sww)),
  c(2025L, 3L),
  info = "us_fiscal_sww ends in 2025 Q3."
)
expect_identical(
  colnames(us_fiscal_sww),
  c("ttr", "gs", "gdp", "FFR", "cons", "rw", "inv", "m2", "ppiic", "pi"),
  info = "us_fiscal_sww has the expected column names."
)
expect_true(
  is.numeric(us_fiscal_sww),
  info = "us_fiscal_sww is numeric."
)
expect_true(
  all(is.finite(us_fiscal_sww)),
  info = "us_fiscal_sww has only finite values."
)

data(us_macro_financial_sw)

expect_true(
  inherits(us_macro_financial_sw, "ts"),
  info = "us_macro_financial_sw is a time series object."
)
expect_identical(
  dim(us_macro_financial_sw),
  c(621L, 10L),
  info = "us_macro_financial_sw has the expected dimensions."
)
expect_equal(
  frequency(us_macro_financial_sw),
  12,
  info = "us_macro_financial_sw is monthly."
)
expect_identical(
  as.integer(start(us_macro_financial_sw)),
  c(1974L, 1L),
  info = "us_macro_financial_sw starts in 1974-01."
)
expect_identical(
  as.integer(end(us_macro_financial_sw)),
  c(2025L, 9L),
  info = "us_macro_financial_sw ends in 2025-09."
)
expect_identical(
  colnames(us_macro_financial_sw),
  c("ip", "p", "R", "hhc", "bc", "GZ", "ES", "TS", "m", "pcm"),
  info = "us_macro_financial_sw has the expected column names."
)
expect_true(
  is.numeric(us_macro_financial_sw),
  info = "us_macro_financial_sw is numeric."
)
expect_true(
  all(is.finite(us_macro_financial_sw)),
  info = "us_macro_financial_sw has only finite values."
)
