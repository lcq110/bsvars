one_horizon = structure(
  list(forecasts = array(c(1, 2, 3), c(1, 1, 3))),
  class = "Forecasts"
)
invisible(capture.output({
  one_horizon_summary = summary(one_horizon)
}))

expect_equal(
  dim(one_horizon_summary[[1]]),
  c(1, 4),
  info = "Forecast summaries retain the horizon dimension when H = 1."
)
expect_equal(
  one_horizon_summary[[1]][1, "mean"],
  2,
  info = "Forecast summaries aggregate posterior draws when H = 1."
)

one_draw = structure(
  list(forecasts = array(c(1, 2, 3), c(1, 3, 1))),
  class = "Forecasts"
)
invisible(capture.output({
  one_draw_summary = summary(one_draw)
}))

expect_equal(
  dim(one_draw_summary[[1]]),
  c(3, 4),
  info = "Forecast summaries retain the horizon dimension when S = 1."
)
expect_equal(
  unname(one_draw_summary[[1]][, "mean"]),
  c(1, 2, 3),
  info = "Forecast summaries retain each horizon when S = 1."
)
