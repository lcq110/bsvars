#include "multiprecision.h"

#include <boost/multiprecision/cpp_bin_float.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

using Backend = boost::multiprecision::cpp_bin_float<
  256,
  boost::multiprecision::backends::digit_base_2
>;
using HighPrecision = boost::multiprecision::number<
  Backend,
  boost::multiprecision::et_off
>;

class PrecisionFactor {
 public:
  explicit PrecisionFactor(const arma::mat& design)
      : rows_(design.n_rows),
        columns_(design.n_cols),
        design_(),
        lower_(columns_ * columns_, HighPrecision(0)) {
    design_.reserve(design.n_elem);
    for (arma::uword column = 0; column < design.n_cols; ++column) {
      for (arma::uword row = 0; row < design.n_rows; ++row) {
        design_.emplace_back(design(row, column));
      }
    }

    std::vector<HighPrecision> precision(
      columns_ * columns_, HighPrecision(0)
    );
    for (std::size_t column = 0; column < columns_; ++column) {
      for (std::size_t other = 0; other <= column; ++other) {
        HighPrecision value = 0;
        for (std::size_t row = 0; row < rows_; ++row) {
          value += at(row, column) * at(row, other);
        }
        precision[column + columns_ * other] = value;
        precision[other + columns_ * column] = value;
      }
    }
    lower_ = cholesky(precision, columns_);
  }

  arma::vec draw(const arma::vec& response) const {
    std::vector<HighPrecision> converted_response;
    converted_response.reserve(response.n_elem);
    for (const double value : response) {
      converted_response.emplace_back(value);
    }
    std::vector<HighPrecision> information(columns_, HighPrecision(0));
    for (std::size_t column = 0; column < columns_; ++column) {
      HighPrecision value = 0;
      for (std::size_t row = 0; row < rows_; ++row) {
        value += at(row, column) * converted_response[row];
      }
      information[column] = value;
    }

    std::vector<HighPrecision> intermediate(columns_, HighPrecision(0));
    forward_solve(lower_, information, intermediate, columns_);

    arma::vec noise(columns_, arma::fill::randn);
    for (std::size_t row = 0; row < columns_; ++row) {
      intermediate[row] += HighPrecision(noise(row));
    }

    std::vector<HighPrecision> coefficients(columns_, HighPrecision(0));
    backward_solve(lower_, intermediate, coefficients, columns_);
    arma::vec result(columns_);
    for (std::size_t row = 0; row < columns_; ++row) {
      result(row) = coefficients[row].convert_to<double>();
    }
    return result;
  }

  arma::mat inverse_cholesky(const int posterior_nu) const {
    std::vector<HighPrecision> inverse_lower(
      columns_ * columns_, HighPrecision(0)
    );
    std::vector<HighPrecision> right_hand_side(
      columns_, HighPrecision(0)
    );
    std::vector<HighPrecision> solution(columns_, HighPrecision(0));
    for (std::size_t column = 0; column < columns_; ++column) {
      std::fill(right_hand_side.begin(), right_hand_side.end(), HighPrecision(0));
      std::fill(solution.begin(), solution.end(), HighPrecision(0));
      right_hand_side[column] = 1;
      forward_solve(lower_, right_hand_side, solution, columns_);
      for (std::size_t row = 0; row < columns_; ++row) {
        inverse_lower[row + columns_ * column] = solution[row];
      }
    }

    std::vector<HighPrecision> covariance(
      columns_ * columns_, HighPrecision(0)
    );
    const HighPrecision scale(posterior_nu);
    for (std::size_t column = 0; column < columns_; ++column) {
      for (std::size_t row = 0; row <= column; ++row) {
        HighPrecision value = 0;
        for (std::size_t k = 0; k < columns_; ++k) {
          value += inverse_lower[k + columns_ * row] *
                   inverse_lower[k + columns_ * column];
        }
        value *= scale;
        covariance[row + columns_ * column] = value;
        covariance[column + columns_ * row] = value;
      }
    }

    const std::vector<HighPrecision> covariance_lower =
      cholesky(covariance, columns_);
    arma::mat result(columns_, columns_, arma::fill::zeros);
    for (std::size_t row = 0; row < columns_; ++row) {
      for (std::size_t column = row; column < columns_; ++column) {
        result(row, column) =
          covariance_lower[column + columns_ * row].convert_to<double>();
      }
    }
    return result;
  }

 private:
  const HighPrecision& at(
      const std::size_t row,
      const std::size_t column) const {
    return design_[row + rows_ * column];
  }

  static std::vector<HighPrecision> cholesky(
      const std::vector<HighPrecision>& matrix,
      const std::size_t dimension) {
    std::vector<HighPrecision> lower(
      dimension * dimension,
      HighPrecision(0)
    );
    for (std::size_t row = 0; row < dimension; ++row) {
      for (std::size_t column = 0; column <= row; ++column) {
        HighPrecision value = matrix[row + dimension * column];
        for (std::size_t k = 0; k < column; ++k) {
          value -= lower[row + dimension * k] *
                   lower[column + dimension * k];
        }
        if (row == column) {
          if (value <= 0) {
            throw std::runtime_error(
              "multiprecision Cholesky is not positive definite"
            );
          }
          lower[row + dimension * column] = sqrt(value);
        } else {
          lower[row + dimension * column] =
            value / lower[column + dimension * column];
        }
      }
    }
    return lower;
  }

  static void forward_solve(
      const std::vector<HighPrecision>& lower,
      const std::vector<HighPrecision>& right_hand_side,
      std::vector<HighPrecision>& solution,
      const std::size_t dimension) {
    for (std::size_t row = 0; row < dimension; ++row) {
      HighPrecision value = right_hand_side[row];
      for (std::size_t column = 0; column < row; ++column) {
        value -= lower[row + dimension * column] * solution[column];
      }
      solution[row] = value / lower[row + dimension * row];
    }
  }

  static void backward_solve(
      const std::vector<HighPrecision>& lower,
      const std::vector<HighPrecision>& right_hand_side,
      std::vector<HighPrecision>& solution,
      const std::size_t dimension) {
    for (std::size_t reverse = 0; reverse < dimension; ++reverse) {
      const std::size_t row = dimension - reverse - 1;
      HighPrecision value = right_hand_side[row];
      for (std::size_t column = row + 1; column < dimension; ++column) {
        value -= lower[column + dimension * row] * solution[column];
      }
      solution[row] = value / lower[row + dimension * row];
    }
  }

  std::size_t rows_;
  std::size_t columns_;
  std::vector<HighPrecision> design_;
  std::vector<HighPrecision> lower_;
};

}  // namespace

namespace bsvars {

arma::vec draw_regression_coefficients_multiprecision(
    const arma::mat& design,
    const arma::vec& response) {
  return PrecisionFactor(design).draw(response);
}

arma::mat chol_inverse_precision_multiprecision(
    const arma::mat& design,
    const int posterior_nu) {
  return PrecisionFactor(design).inverse_cholesky(posterior_nu);
}

}  // namespace bsvars
