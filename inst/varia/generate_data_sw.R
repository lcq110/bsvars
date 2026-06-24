# Provenance script for data/us_macro_financial_sw.rda.
#
# The packaged object was generated from the author-prepared processed file:
#   /Users/felix/Documents/GitHub/bsvar_hmsh/BPSSdata/Forecasting study/data/finaldata_processed.rda
#
# That file is generated from:
#   /Users/felix/Documents/GitHub/bsvar_hmsh/empirical/finaldata.rda
#
# which is equivalent to:
#   /Users/felix/Documents/GitHub/bsvar_hmsh/BPSSdata/data/finaldata.rda
#
# The upstream recipe is:
#   /Users/felix/Documents/GitHub/bsvar_hmsh/BPSSdata/data/finaldata_process.R
#
# This script records the final transformation from that local source object to
# the package data object. It is kept under inst/varia, which is excluded from
# the CRAN build by .Rbuildignore.

source_path <- "/Users/felix/Documents/GitHub/bsvar_hmsh/BPSSdata/Forecasting study/data/finaldata_processed.rda"
if (!file.exists(source_path)) {
  stop("Cannot find the author-prepared finaldata_processed.rda source file.")
}

load(source_path) # loads df and stationary
stopifnot(exists("df"))

required_columns <- c("DATE", "ip", "p", "r", "hhc", "bc", "gz", "es", "ts", "m1", "pcm")
stopifnot(all(required_columns %in% colnames(df)))

df <- df[order(df$DATE), ]
dates <- as.Date(df$DATE)
sample_rows <- dates >= as.Date("1974-01-01") & dates <= as.Date("2025-09-01")

transformed <- cbind(
  ip = df$ip,
  p = df$p,
  R = df$r,
  hhc = df$hhc,
  bc = df$bc,
  GZ = df$gz,
  ES = df$es,
  TS = df$ts,
  m = df$m1,
  pcm = df$pcm
)

us_macro_financial_sw <- ts(
  transformed[sample_rows, ],
  start = c(1974, 1),
  frequency = 12
)

stopifnot(identical(dim(us_macro_financial_sw), c(621L, 10L)))
stopifnot(identical(colnames(us_macro_financial_sw), c("ip", "p", "R", "hhc", "bc", "GZ", "ES", "TS", "m", "pcm")))
stopifnot(all(is.finite(us_macro_financial_sw)))

save(us_macro_financial_sw, file = "data/us_macro_financial_sw.rda")
