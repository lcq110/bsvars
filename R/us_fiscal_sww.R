#' @title A 10-variable US fiscal and macroeconomic forecasting system for the period 1959 Q2 -- 2025 Q3
#'
#' @description A quarterly US fiscal and macroeconomic system used in the
#' International Journal of Forecasting paper by Fei Shang, Xiaolei (Adam) Wang,
#' and Tomasz Woźniak.
#'
#' @usage data(us_fiscal_sww)
#'
#' @format A matrix and a \code{ts} object with 266 quarterly observations on
#' 10 variables:
#' \describe{
#'   \item{ttr}{US total tax revenue, real per person growth rate}
#'   \item{gs}{US government spending, real per person growth rate}
#'   \item{gdp}{US gross domestic product, real per person growth rate}
#'   \item{FFR}{Federal Funds Effective Rate}
#'   \item{cons}{private consumption, real per person growth rate}
#'   \item{rw}{real compensation per hour growth rate}
#'   \item{inv}{private non-residential investment, real per person growth rate}
#'   \item{m2}{M2 monetary aggregate growth rate}
#'   \item{ppiic}{producer price index for industrial commodities growth rate}
#'   \item{pi}{GDP deflator inflation measure}
#' }
#'
#' The first three fiscal variables follow the construction in Mertens and Ravn
#' (2014). The 10-variable system extends this fiscal block with the variables
#' used by Mountford and Uhlig (2009).
#'
#' @references
#' Mountford, A., and Uhlig, H. (2009) What are the effects of fiscal policy
#' shocks? \emph{Journal of Applied Econometrics}, 24(6), 960--992,
#' DOI: \doi{10.1002/jae.1079}.
#'
#' Mertens, K., and Ravn, M.O. (2014) A Reconciliation of SVAR and Narrative
#' Estimates of Tax Multipliers, \emph{Journal of Monetary Economics}, 68(S),
#' S1--S19. DOI: \doi{10.1016/j.jmoneco.2013.04.004}.
#'
#' @source
#' U.S. Bureau of Economic Analysis, National Income and Product Accounts,
#' \url{https://www.bea.gov/}
#'
#' FRED Economic Database, Federal Reserve Bank of St. Louis,
#' \url{https://fred.stlouisfed.org/}
#'
#' @examples
#' data(us_fiscal_sww)   # upload the data
#' plot(us_fiscal_sww)   # plot the data
"us_fiscal_sww"
