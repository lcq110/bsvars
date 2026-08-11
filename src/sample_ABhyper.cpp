

#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include "progress.hpp"
#include "Rcpp/Rmath.h"

#include "utils.h"
#include "multiprecision.h"

using namespace Rcpp;
using namespace arma;


namespace {

arma::vec draw_regression_coefficients (
    const arma::mat&  prior_chol,
    const arma::vec&  prior_mean,
    const arma::mat&  weighted_observations,
    const arma::vec&  weighted_response,
    const arma::mat&  restrictions,
    const double      prior_scale
) {
  mat design = join_cols(
    prior_scale * prior_chol * trans(restrictions),
    weighted_observations * trans(restrictions)
  );
  vec response = join_cols(
    prior_scale * prior_chol * prior_mean,
    weighted_response
  );
  return bsvars::draw_regression_coefficients_multiprecision(
    design,
    response
  );
}


arma::mat chol_inverse_precision (
    const arma::mat&  prior_chol,
    const arma::mat&  weighted_observations,
    const arma::mat&  restrictions,
    const double      prior_scale,
    const int         posterior_nu
) {
  const int dimension = restrictions.n_rows;
  mat covariance_chol;
  mat design = join_cols(
    prior_scale * prior_chol * trans(restrictions),
    trans(weighted_observations) * trans(restrictions)
  );
  mat Q;
  mat precision_chol;
  qr_econ(Q, precision_chol, design);
  precision_chol = trimatu(precision_chol);

  const double reciprocal_condition = rcond(precision_chol);
  const double threshold =
    dimension * std::sqrt(std::numeric_limits<double>::epsilon());
  if (!(reciprocal_condition > threshold)) {
    return bsvars::chol_inverse_precision_multiprecision(
      design,
      posterior_nu
    );
  }

  mat inverse_factor = solve(
    trimatl(trans(precision_chol)),
    eye<mat>(dimension, dimension),
    solve_opts::no_approx
  );
  qr_econ(Q, covariance_chol, inverse_factor);
  covariance_chol = trimatu(covariance_chol);
  for (int i=0; i<dimension; i++) {
    if (covariance_chol(i,i) < 0) {
      covariance_chol.row(i) *= -1;
    }
  }

  return sqrt(posterior_nu) * covariance_chol;
}

} // anonymous namespace



/*______________________function sample_A_homosk1______________________*/
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_A_homosk1 (
    arma::mat&        aux_A,          // NxK
    const arma::mat&  aux_B,          // NxN
    const arma::mat&  aux_hyper,      // (2*N+1) x 2 :: col 0 for B, col 1 for A
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VA  // restrictions on A
) {
  // the function changes the value of aux_A by reference
  const int N         = aux_A.n_rows;
  const int K         = aux_A.n_cols;
  
  mat prior_A_mean    = as<mat>(prior["A"]);
  mat prior_A_Vinv    = as<mat>(prior["A_V_inv"]);
  mat prior_A_chol    = trimatu(chol(prior_A_Vinv));
  rowvec    zerosA(K);
  
  for (int n=0; n<N; n++) {
    mat   A0          = aux_A;
    A0.row(n)         = zerosA;
    vec   zn          = vectorise( aux_B * (Y - A0 * X) );
    mat   Wn          = kron( trans(X), aux_B.col(n) );
    
    vec draw = draw_regression_coefficients(
      prior_A_chol,
      trans(prior_A_mean.row(n)),
      Wn,
      zn,
      VA(n),
      pow(aux_hyper(n,1), -0.5)
    );
    aux_A.row(n)      = trans(draw) * VA(n);
  } // END n loop
  
  return aux_A;
} // END sample_A_homosk1



/*______________________function sample_A_heterosk1 ______________________*/
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_A_heterosk1 (
    arma::mat&        aux_A,          // NxK
    const arma::mat&  aux_B,          // NxN
    const arma::mat&  aux_hyper,      // (2*N+1) x 2 :: col 0 for B, col 1 for A
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VA  // restrictions on A
) {
  // the function changes the value of aux_A by reference
  const int N         = aux_A.n_rows;
  const int K         = aux_A.n_cols;
  
  mat prior_A_mean    = as<mat>(prior["A"]);
  mat prior_A_Vinv    = as<mat>(prior["A_V_inv"]);
  mat prior_A_chol    = trimatu(chol(prior_A_Vinv));
  rowvec    zerosA(K);
  vec sigma_vectorised= vectorise(aux_sigma);
  
  for (int n=0; n<N; n++) {
    mat   A0          = aux_A;
    A0.row(n)         = zerosA;
    vec   zn          = vectorise( aux_B * (Y - A0 * X) );
    mat   zn_sigma    = zn / sigma_vectorised;
    mat   Wn          = kron( trans(X), aux_B.col(n) );
    mat   Wn_sigma    = Wn.each_col() / sigma_vectorised;
    
    vec draw = draw_regression_coefficients(
      prior_A_chol,
      trans(prior_A_mean.row(n)),
      Wn_sigma,
      zn_sigma,
      VA(n),
      pow(aux_hyper(n,1), -0.5)
    );
    aux_A.row(n)      = trans(draw) * VA(n);
  } // END n loop
  
  return aux_A;
} // END sample_A_heterosk1



/*______________________function sample_B_homosk1______________________*/
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_B_homosk1 (
    arma::mat&        aux_B,          // NxN
    const arma::mat&  aux_A,          // NxK
    const arma::mat&  aux_hyper,      // (2*N+1) x 2 :: col 0 for B, col 1 for A
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
) {
  // the function changes the value of aux_B by reference
  const int N               = aux_B.n_rows;
  const int T               = Y.n_cols;
  
  const int posterior_nu    = T + as<int>(prior["B_nu"]);
  // The generalized-normal determinant power is posterior_nu - N.
  const int radial_degrees  = posterior_nu - N + 1;
  mat prior_SS_inv          = as<mat>(prior["B_V_inv"]);
  mat prior_SS_chol         = trimatu(chol(prior_SS_inv));
  mat shocks                = Y - aux_A * X;
  
  for (int n=0; n<N; n++) {
    // sample B
    mat Un                  = chol_inverse_precision(
      prior_SS_chol,
      shocks,
      VB(n),
      pow(aux_hyper(n,0), -0.5),
      posterior_nu
    );
    mat B_tmp               = aux_B;
    B_tmp.shed_row(n);
    rowvec w                = trans(orthogonal_complement_matrix_TW(B_tmp.t()));
    vec w1_tmp              = trans(w * VB(n).t() * Un.t());
    double w1w1_tmp         = as_scalar(sum(pow(w1_tmp, 2)));
    mat w1                  = w1_tmp.t()/sqrt(w1w1_tmp);
    mat Wn;
    const int rn            = VB(n).n_rows;
    if (rn==1) {
      Wn                    = w1;
    } else {
      Wn                    = join_rows(w1.t(), orthogonal_complement_matrix_TW(w1.t()));
    }
    
    vec   alpha(rn);
    vec   u(radial_degrees, fill::randn);
    u                      *= pow(posterior_nu, -0.5);
    alpha(0)                = sqrt(as_scalar(sum(square(u))));
    if (R::runif(0,1)<0.5) {
      alpha(0)       *= -1;
    }
    if (rn>1){
      vec nn(rn-1, fill::randn);
      nn                   *= pow(posterior_nu, -0.5);
      alpha.rows(1,rn-1)    = nn;
    }
    rowvec b0n              = alpha.t() * Wn * Un;
    aux_B.row(n)            = b0n * VB(n);
  } // END n loop
  
  return aux_B;
} // END sample_B_homosk1



/*______________________function sample_B_heterosk1______________________*/
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_B_heterosk1 (
    arma::mat&        aux_B,          // NxN
    const arma::mat&  aux_A,          // NxK
    const arma::mat&  aux_hyper,      // (2*N+1) x 2 :: col 0 for B, col 1 for A
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
) {
  // the function changes the value of aux_B0 and aux_Bplus by reference (filling it with a new draw)
  const int N               = aux_B.n_rows;
  const int T               = Y.n_cols;
  
  const int posterior_nu    = T + as<int>(prior["B_nu"]);
  // The generalized-normal determinant power is posterior_nu - N.
  const int radial_degrees  = posterior_nu - N + 1;
  mat prior_SS_inv          = as<mat>(prior["B_V_inv"]);
  mat prior_SS_chol         = trimatu(chol(prior_SS_inv));
  mat shocks                = Y - aux_A * X;
  
  
  for (int n=0; n<N; n++) {
    
    // set scale matrix
    mat shocks_sigma        = shocks.each_row() / aux_sigma.row(n);
    
    // sample B
    mat Un                  = chol_inverse_precision(
      prior_SS_chol,
      shocks_sigma,
      VB(n),
      pow(aux_hyper(n,0), -0.5),
      posterior_nu
    );
    mat B_tmp               = aux_B;
    B_tmp.shed_row(n);
    rowvec w                = trans(orthogonal_complement_matrix_TW(B_tmp.t()));
    vec w1_tmp              = trans(w * VB(n).t() * Un.t());
    double w1w1_tmp         = as_scalar(sum(pow(w1_tmp, 2)));
    mat w1                  = w1_tmp.t()/sqrt(w1w1_tmp);
    mat Wn;
    const int rn            = VB(n).n_rows;
    if (rn==1) {
      Wn                    = w1;
    } else {
      Wn                    = join_rows(w1.t(), orthogonal_complement_matrix_TW(w1.t()));
    }
    
    vec   alpha(rn);
    vec   u(radial_degrees, fill::randn);
    u                      *= pow(posterior_nu, -0.5);
    alpha(0)                = pow(as_scalar(sum(pow(u,2))), 0.5);
    if (R::runif(0,1)<0.5) {
      alpha(0)       *= -1;
    }
    if (rn>1){
      vec nn(rn-1, fill::randn);
      nn                   *= pow(posterior_nu, -0.5);
      alpha.rows(1,rn-1)    = nn;
    }
    rowvec b0n              = alpha.t() * Wn * Un;
    aux_B.row(n)           = b0n * VB(n);
  } // END n loop
  
  return aux_B;
} // END sample_B_heterosk1



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_hyperparameters (
    arma::mat&              aux_hyper,       // (2*N+1) x 2 :: col 0 for B, col 1 for A
    const arma::mat&        aux_B,            // NxN
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::field<arma::mat>& VA,
    const Rcpp::List&       prior
) {
  // the function returns aux_hyper by reference (filling it with a new draw)
  
  const int N = aux_B.n_rows;
  
  double prior_B_nu           = as<double>(prior["B_nu"]);
  double prior_hyper_nu_B     = as<double>(prior["hyper_nu_B"]);
  double prior_hyper_a_B      = as<double>(prior["hyper_a_B"]);
  double prior_hyper_s_BB     = as<double>(prior["hyper_s_BB"]);
  double prior_hyper_nu_BB    = as<double>(prior["hyper_nu_BB"]);
  
  double prior_hyper_nu_A     = as<double>(prior["hyper_nu_A"]);
  double prior_hyper_a_A      = as<double>(prior["hyper_a_A"]);
  double prior_hyper_s_AA     = as<double>(prior["hyper_s_AA"]);
  double prior_hyper_nu_AA    = as<double>(prior["hyper_nu_AA"]);
  
  mat   prior_A               = as<mat>(prior["A"]);
  mat   prior_A_V_inv         = as<mat>(prior["A_V_inv"]);
  mat   prior_B_V_inv         = as<mat>(prior["B_V_inv"]);
  
  // aux_B - related hyper-parameters 
  vec     ss_tmp      = aux_hyper.submat(N, 0, 2 * N - 1, 0);
  double  scale_tmp   = prior_hyper_s_BB + 2 * sum(ss_tmp);
  double  shape_tmp   = prior_hyper_nu_BB + 2 * N * prior_hyper_a_B;
  aux_hyper(2 * N, 0) = scale_tmp / R::rchisq(shape_tmp);
  
  // aux_A - related hyper-parameters 
  ss_tmp              = aux_hyper.submat(N, 1, 2 * N - 1, 1);
  scale_tmp           = prior_hyper_s_AA + 2 * sum(ss_tmp);
  shape_tmp           = prior_hyper_nu_AA + 2 * N * prior_hyper_a_A;
  aux_hyper(2 * N, 1) = scale_tmp / R::rchisq(shape_tmp);
  
  for (int n=0; n<N; n++) {
    
    // count unrestricted elements of aux_B's row
    int rn            = VB(n).n_rows;
    int rna           = VA(n).n_rows;
    
    // aux_B - related hyper-parameters 
    scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 0))) + (1 / aux_hyper(2 * N, 0)));
    shape_tmp         = prior_hyper_a_B + prior_hyper_nu_B / 2;
    aux_hyper(N + n, 0) = R::rgamma(shape_tmp, scale_tmp);
    
    scale_tmp         = aux_hyper(N + n, 0) + as_scalar(aux_B.row(n) * prior_B_V_inv * aux_B.row(n).t());
    shape_tmp         = prior_hyper_nu_B + rn + prior_B_nu - N;
    aux_hyper(n, 0)   = scale_tmp / R::rchisq(shape_tmp);
    
    // aux_A - related hyper-parameters 
    scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 1))) + (1 / aux_hyper(2 * N, 1)));
    shape_tmp         = prior_hyper_a_A + prior_hyper_nu_A / 2;
    aux_hyper(N + n, 1) = R::rgamma(shape_tmp, scale_tmp);
    
    scale_tmp         = aux_hyper(N + n, 1) + 
      as_scalar((aux_A.row(n) - prior_A.row(n)) * prior_A_V_inv * trans(aux_A.row(n) - prior_A.row(n)));
    shape_tmp         = prior_hyper_nu_A + rna;
    aux_hyper(n, 1)   = scale_tmp / R::rchisq(shape_tmp);
  } // END n loop
  
  return aux_hyper;
} // END sample_hyperparameters
