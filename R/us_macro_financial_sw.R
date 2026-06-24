#' @title A 10-variable US monthly macro-financial system for the period 1974-01 -- 2025-09
#'
#' @description A monthly US macro-financial system used in Shang and Woźniak
#' (2026). The data contain the variables used in the empirical application on
#' sparse heterogeneous Markov switching heteroskedasticity.
#'
#' @usage data(us_macro_financial_sw)
#'
#' @format A matrix and a \code{ts} object with 621 monthly observations on
#' 10 variables:
#' \describe{
#'   \item{ip}{industrial production, 12-month log difference}
#'   \item{p}{PCE price index, 12-month log difference}
#'   \item{R}{Federal Funds Effective Rate, in decimal form}
#'   \item{hhc}{household credit, 12-month log difference}
#'   \item{bc}{business credit, 12-month log difference}
#'   \item{GZ}{corporate bond spread, in decimal form}
#'   \item{ES}{external finance spread, in decimal form}
#'   \item{TS}{term spread, in decimal form}
#'   \item{m}{M2 monetary aggregate, 12-month log difference}
#'   \item{pcm}{Bloomberg Commodity Index, 12-month log difference}
#' }
#'
#' Industrial production, prices, household credit, business credit, money, and
#' commodity prices are expressed as 12-month log differences. Interest rates and
#' spreads are expressed in decimal form.
#'
#' @references
#' Shang, F., and Woźniak, T. (2026) Identification Verification for Structural
#' Vector Autoregressions with Sparse Heterogeneous Markov Switching
#' Heteroskedasticity. arXiv:2603.16035, \doi{10.48550/arXiv.2603.16035}.
#'
#' @source
#' Author-prepared monthly macro-financial dataset from the Shang and Woźniak
#' (2026) empirical application, based on FRED macro-financial series, corporate
#' bond spread data, spliced commercial paper spreads, and the Bloomberg
#' Commodity Index.
#'
#' @examples
#' data(us_macro_financial_sw)   # upload the data
#' plot(us_macro_financial_sw)   # plot the data
"us_macro_financial_sw"
