
#ifndef _UTILS_H_
#define _UTILS_H_

#include "bsvars_RcppArmadillo.h"

arma::mat orthogonal_complement_matrix_TW (const arma::mat& x);

arma::vec log_mean (arma::mat log_density);

std::string ordinal(int n);

#endif  // _UTILS_H_