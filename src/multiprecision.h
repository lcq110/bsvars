#ifndef BSVARS_MULTIPRECISION_H
#define BSVARS_MULTIPRECISION_H

#include <RcppArmadillo.h>

namespace bsvars {

arma::vec draw_regression_coefficients_multiprecision(
  const arma::mat& design,
  const arma::vec& response
);

arma::mat chol_inverse_precision_multiprecision(
  const arma::mat& design,
  int posterior_nu
);

}  // namespace bsvars

#endif
