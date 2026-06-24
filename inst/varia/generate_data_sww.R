# Provenance script for data/us_fiscal_sww.rda.
#
# The packaged object was copied unchanged from:
#   /Users/felix/Documents/GitHub/bsvarsFORE/data/us_fiscal_sww.rda
#
# Original recipe:
#   /Users/felix/Documents/GitHub/bsvarsFORE/data/us_fiscal_sww.R
#
# This script records the data construction used in that source project. It is
# kept under inst/varia, which is excluded from the CRAN build by .Rbuildignore.
# It is not a runtime dependency of the package.

stopifnot(requireNamespace("fredr", quietly = TRUE))
stopifnot(requireNamespace("xts", quietly = TRUE))
stopifnot(requireNamespace("zoo", quietly = TRUE))

path <- "inst/varia/"

# Population, converted from monthly FRED data to quarterly averages.
pop_tmp <- fredr::fredr("CNP16OV")
pop <- xts::to.quarterly(xts::xts(pop_tmp$value, pop_tmp$date), OHLC = FALSE)
pop <- log(pop)

# GDP deflator: BEA NIPA Table 1.1.9, line 1.
pi_read <- read.csv(file.path(path, "Table1.1.9.csv"), header = FALSE)
pi_value <- log(as.numeric(t(pi_read[6, 3:ncol(pi_read)])))
pi_date <- zoo::as.yearqtr(paste(t(pi_read[4, 3:ncol(pi_read)]), t(pi_read[5, 3:ncol(pi_read)])))
pi <- xts::xts(pi_value, pi_date)

# GDP: BEA NIPA Table 1.1.5, line 1.
gdp_read <- read.csv(file.path(path, "Table1.1.5.csv"), header = FALSE)
gdp_value <- log(as.numeric(t(gdp_read[6, 3:ncol(gdp_read)])))
gdp_date <- zoo::as.yearqtr(paste(t(gdp_read[4, 3:ncol(gdp_read)]), t(gdp_read[5, 3:ncol(gdp_read)])))
gdp <- xts::xts(gdp_value, gdp_date)
gdp <- gdp - pop - pi
gdp <- 100 * na.omit(diff(gdp))

# Government spending: BEA NIPA Table 3.9.5, line 9 in the current layout.
gs_read <- read.csv(file.path(path, "Table3.9.5.csv"), header = FALSE)
gs_value <- log(as.numeric(t(gs_read[14, 3:ncol(gs_read)])))
gs_date <- zoo::as.yearqtr(paste(t(gs_read[4, 3:ncol(gs_read)]), t(gs_read[5, 3:ncol(gs_read)])))
gs <- xts::xts(gs_value, gs_date)
gs <- gs - pop - pi
gs <- 100 * na.omit(diff(gs))

# Total tax revenue:
# Federal current tax receipts plus contributions for government social
# insurance, less corporate income taxes from Federal Reserve Banks.
ttr_read <- read.csv(file.path(path, "Table3.2.csv"), header = FALSE)
ttr_ctr <- as.numeric(t(ttr_read[7, 3:ncol(ttr_read)]))
ttr_cgsi <- as.numeric(t(ttr_read[15, 3:ncol(ttr_read)]))
ttr_cit <- as.numeric(t(ttr_read[13, 3:ncol(ttr_read)]))
ttr_value <- log(ttr_ctr + ttr_cgsi - ttr_cit)
ttr_date <- zoo::as.yearqtr(paste(t(ttr_read[4, 3:ncol(ttr_read)]), t(ttr_read[5, 3:ncol(ttr_read)])))
ttr <- xts::xts(ttr_value, ttr_date)
ttr <- ttr - pop - pi
ttr <- 100 * na.omit(diff(ttr))

us_fiscal_sml <- na.omit(cbind(ttr, gs, gdp))
colnames(us_fiscal_sml) <- c("ttr", "gs", "gdp")

# Mountford and Uhlig (2009) extension variables.
cons_read <- read.csv(file.path(path, "Table1.1.5.csv"), header = FALSE)
cons_value <- log(as.numeric(t(cons_read[7, 3:ncol(cons_read)])))
cons_date <- zoo::as.yearqtr(paste(t(cons_read[4, 3:ncol(cons_read)]), t(cons_read[5, 3:ncol(cons_read)])))
cons <- xts::xts(cons_value, cons_date)
cons <- cons - pop - pi
cons <- 100 * na.omit(diff(cons))
colnames(cons) <- "cons"

rw_tmp <- fredr::fredr("COMPRNFB")
rw <- xts::xts(log(rw_tmp$value), zoo::as.yearqtr(rw_tmp$date))
rw <- 100 * na.omit(diff(rw))
colnames(rw) <- "rw"

pnri_read <- read.csv(file.path(path, "Table1.1.5.csv"), header = FALSE)
gpdi_value <- as.numeric(t(pnri_read[12, 3:ncol(pnri_read)]))
ri_value <- as.numeric(t(pnri_read[18, 3:ncol(pnri_read)]))
inv_date <- zoo::as.yearqtr(paste(t(pnri_read[4, 3:ncol(pnri_read)]), t(pnri_read[5, 3:ncol(pnri_read)])))
inv <- xts::xts(log(gpdi_value - ri_value), inv_date)
inv <- inv - pop - pi
inv <- 100 * na.omit(diff(inv))
colnames(inv) <- "inv"

ffr_tmp <- fredr::fredr("FEDFUNDS")
FFR <- xts::to.quarterly(xts::xts(ffr_tmp$value, ffr_tmp$date), OHLC = FALSE)
FFR <- xts::xts(FFR, zoo::as.yearqtr(zoo::index(FFR)))
colnames(FFR) <- "FFR"

m_tmp <- fredr::fredr("M2SL")
m <- xts::to.quarterly(xts::xts(m_tmp$value, m_tmp$date), OHLC = FALSE)
m <- 100 * na.omit(diff(log(m)))
colnames(m) <- "m2"

ppiic_tmp <- fredr::fredr("PPIIDC")
ppiic <- xts::to.quarterly(xts::xts(ppiic_tmp$value, ppiic_tmp$date), OHLC = FALSE)
ppiic <- 100 * na.omit(diff(log(ppiic)))
colnames(ppiic) <- "ppiic"

pi <- 100 * na.omit(diff(log(pi)))
colnames(pi) <- "pi"

us_fiscal_sww <- ts(
  na.omit(cbind(us_fiscal_sml, FFR, cons, rw, inv, m, ppiic, pi)),
  start = c(1959, 2),
  frequency = 4
)

save(us_fiscal_sww, file = "data/us_fiscal_sww.rda")
