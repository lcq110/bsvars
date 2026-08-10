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

one_draw$forecast_factor = list(
  representation = "B_sigma2_v1",
  B = array(1, c(1, 1, 1)),
  sigma2 = array(1, c(1, 3, 1))
)
invisible(capture.output({
  factor_summary = summary(one_draw)
}))

expect_identical(
  factor_summary,
  one_draw_summary,
  info = "Forecast summaries ignore the additional structural-factor field."
)

one_draw$Y = matrix(1:4, 1, 4)
plot_result = local({
  plot_file = tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    grDevices::dev.off()
    unlink(plot_file)
  })
  plot(one_draw)
})

expect_identical(
  plot_result,
  one_draw,
  info = "Forecast plots ignore the additional structural-factor field."
)
