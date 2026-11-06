set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvar$new(us_fiscal_lsuw)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)

for (standardise in c(FALSE, TRUE)) {
  irf               <- compute_impulse_responses(run_no1, horizon = 2, standardise = standardise)
  cirf              <- irf
  attr(cirf, "cumulative") <- TRUE
  cirf[,,2,]        <- irf[,,1,] + irf[,,2,]
  cirf[,,3,]        <- irf[,,1,] + irf[,,2,] + irf[,,3,]

  expect_equal(
    compute_cumulative_impulse_responses(run_no1, horizon = 2, standardise = standardise), cirf,
    info = "compute_cumulative_impulse_responses: cumulative sums of impulse responses per draw."
  )
}
