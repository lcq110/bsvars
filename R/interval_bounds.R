# Rows are quantities, columns are posterior draws.
.interval_bounds = function(draws, probability, interval) {
  stopifnot("Argument probability must be a number from interval (0,1)." =
              length(probability) == 1 && probability > 0 && probability < 1)
  if (interval == "hpd") {
    return(coda::HPDinterval(coda::mcmc(t(draws)), prob = probability))
  }
  probs = c(50 - 50 * probability, 50 + 50 * probability) / 100
  t(apply(draws, 1, stats::quantile, probs = probs))
}

.interval_names = function(probability, interval) {
  if (interval == "hpd") {
    return(paste0(100 * probability, "% HPD ", c("lower", "upper")))
  }
  paste0(format(c(50 - 50 * probability, 50 + 50 * probability), trim = TRUE),
         "% quantile")
}
