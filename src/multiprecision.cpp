#include "multiprecision.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
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
using ExactProduct = boost::multiprecision::uint128_t;

struct Binary64 {
  std::uint64_t significand = 0;
  int exponent = 0;
  bool negative = false;
};

struct BinaryRange {
  bool has_nonzero = false;
  int minimum_exponent = 0;
  int maximum_top_exponent = 0;
};

struct ExactDotPlan {
  int base_exponent = 0;
  std::size_t limb_count = 1;
};

Binary64 decompose_binary64(const double value) {
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  const std::uint64_t fraction = bits & 0x000fffffffffffffULL;
  const int encoded_exponent = static_cast<int>((bits >> 52) & 0x7ff);
  if (encoded_exponent == 0x7ff) {
    throw std::runtime_error("non-finite exact dot input");
  }

  Binary64 result;
  result.negative = (bits >> 63) != 0;
  if (encoded_exponent == 0) {
    result.significand = fraction;
    result.exponent = -1074;
  } else {
    result.significand = (std::uint64_t(1) << 52) | fraction;
    result.exponent = encoded_exponent - 1023 - 52;
  }
  return result;
}

std::size_t bit_width(std::uint64_t value) {
  std::size_t width = 0;
  while (value != 0) {
    ++width;
    value >>= 1;
  }
  return width;
}

std::size_t ceiling_log2(const std::size_t value) {
  std::size_t width = 0;
  std::size_t remainder = value > 0 ? value - 1 : 0;
  while (remainder != 0) {
    ++width;
    remainder >>= 1;
  }
  return width;
}

BinaryRange binary_range(const std::vector<Binary64>& values) {
  BinaryRange result;
  result.minimum_exponent = std::numeric_limits<int>::max();
  result.maximum_top_exponent = std::numeric_limits<int>::min();
  for (const Binary64& value : values) {
    if (value.significand == 0) continue;
    result.has_nonzero = true;
    result.minimum_exponent = std::min(
      result.minimum_exponent,
      value.exponent
    );
    result.maximum_top_exponent = std::max(
      result.maximum_top_exponent,
      value.exponent + static_cast<int>(bit_width(value.significand))
    );
  }
  return result;
}

ExactDotPlan make_exact_dot_plan(
    const BinaryRange& left,
    const BinaryRange& right,
    const std::size_t terms) {
  if (!left.has_nonzero || !right.has_nonzero) return ExactDotPlan{};

  ExactDotPlan result;
  result.base_exponent = left.minimum_exponent + right.minimum_exponent;
  const int product_top = left.maximum_top_exponent +
    right.maximum_top_exponent - result.base_exponent;
  const std::size_t required_signed_bits =
    static_cast<std::size_t>(product_top) + ceiling_log2(terms) + 1;
  result.limb_count = (required_signed_bits + 63) / 64;
  return result;
}

class ExactAccumulator {
 public:
  explicit ExactAccumulator(const std::size_t limb_count)
      : limbs_(limb_count, 0) {}

  void add_product(
      const Binary64& left,
      const Binary64& right,
      const int base_exponent) {
    if (left.significand == 0 || right.significand == 0) return;

    const int shift = left.exponent + right.exponent - base_exponent;
    const ExactProduct product = ExactProduct(left.significand) *
      ExactProduct(right.significand);
    const std::uint64_t low = static_cast<std::uint64_t>(product);
    const std::uint64_t high = static_cast<std::uint64_t>(product >> 64);
    const std::size_t first_limb = static_cast<std::size_t>(shift / 64);
    const unsigned offset = static_cast<unsigned>(shift % 64);

    std::array<std::uint64_t, 3> pieces{};
    if (offset == 0) {
      pieces[0] = low;
      pieces[1] = high;
    } else {
      pieces[0] = low << offset;
      pieces[1] = (high << offset) | (low >> (64 - offset));
      pieces[2] = high >> (64 - offset);
    }

    if (left.negative == right.negative) {
      add(first_limb, pieces);
    } else {
      subtract(first_limb, pieces);
    }
  }

  HighPrecision to_high_precision(const int base_exponent) const {
    std::vector<std::uint64_t> magnitude = limbs_;
    const bool negative = (magnitude.back() >> 63) != 0;
    if (negative) {
      std::uint64_t carry = 1;
      for (std::uint64_t& limb : magnitude) {
        const std::uint64_t inverted = ~limb;
        limb = inverted + carry;
        carry = limb < inverted;
      }
    }

    boost::multiprecision::cpp_int exact;
    boost::multiprecision::import_bits(
      exact,
      magnitude.begin(),
      magnitude.end(),
      64,
      false
    );
    if (exact == 0) return HighPrecision(0);
    const unsigned trailing = boost::multiprecision::lsb(exact);
    exact >>= trailing;
    HighPrecision result(exact);
    result = ldexp(result, base_exponent + static_cast<int>(trailing));
    return negative ? -result : result;
  }

 private:
  // Signed headroom makes any final carry or borrow two's-complement wrap.
  void add(
      const std::size_t first,
      const std::array<std::uint64_t, 3>& pieces) {
    std::uint64_t carry = 0;
    std::size_t index = first;
    for (const std::uint64_t piece : pieces) {
      if (index == limbs_.size()) {
        if (piece != 0) {
          throw std::runtime_error("exact accumulator overflow");
        }
        continue;
      }
      const std::uint64_t first_sum = limbs_[index] + piece;
      const std::uint64_t first_carry = first_sum < limbs_[index];
      const std::uint64_t sum = first_sum + carry;
      carry = first_carry | (sum < first_sum);
      limbs_[index] = sum;
      ++index;
    }
    while (carry != 0 && index < limbs_.size()) {
      const std::uint64_t sum = limbs_[index] + carry;
      carry = sum < limbs_[index];
      limbs_[index] = sum;
      ++index;
    }
  }

  void subtract(
      const std::size_t first,
      const std::array<std::uint64_t, 3>& pieces) {
    std::uint64_t borrow = 0;
    std::size_t index = first;
    for (const std::uint64_t piece : pieces) {
      if (index == limbs_.size()) {
        if (piece != 0) {
          throw std::runtime_error("exact accumulator overflow");
        }
        continue;
      }
      const std::uint64_t required = piece + borrow;
      const std::uint64_t required_overflow = required < piece;
      const std::uint64_t current = limbs_[index];
      limbs_[index] = current - required;
      borrow = required_overflow | (current < required);
      ++index;
    }
    while (borrow != 0 && index < limbs_.size()) {
      const std::uint64_t current = limbs_[index];
      limbs_[index] = current - 1;
      borrow = current == 0;
      ++index;
    }
  }

  std::vector<std::uint64_t> limbs_;
};

HighPrecision exact_dot(
    const Binary64* left,
    const Binary64* right,
    const std::size_t size,
    const ExactDotPlan& plan) {
  ExactAccumulator accumulator(plan.limb_count);
  for (std::size_t index = 0; index < size; ++index) {
    accumulator.add_product(left[index], right[index], plan.base_exponent);
  }
  return accumulator.to_high_precision(plan.base_exponent);
}

class PrecisionFactor {
 public:
  explicit PrecisionFactor(const arma::mat& design)
      : rows_(design.n_rows),
        columns_(design.n_cols),
        design_(),
        design_range_(),
        lower_(columns_ * columns_, HighPrecision(0)) {
    design_.reserve(design.n_elem);
    for (arma::uword column = 0; column < design.n_cols; ++column) {
      for (arma::uword row = 0; row < design.n_rows; ++row) {
        design_.push_back(decompose_binary64(design(row, column)));
      }
    }
    design_range_ = binary_range(design_);
    const ExactDotPlan precision_plan = make_exact_dot_plan(
      design_range_,
      design_range_,
      rows_
    );

    std::vector<HighPrecision> precision(
      columns_ * columns_, HighPrecision(0)
    );
    for (std::size_t column = 0; column < columns_; ++column) {
      for (std::size_t other = 0; other <= column; ++other) {
        const HighPrecision value = exact_dot(
          design_.data() + rows_ * column,
          design_.data() + rows_ * other,
          rows_,
          precision_plan
        );
        precision[column + columns_ * other] = value;
        precision[other + columns_ * column] = value;
      }
    }
    lower_ = cholesky(precision, columns_);
  }

  arma::vec draw(const arma::vec& response) const {
    std::vector<Binary64> decomposed_response;
    decomposed_response.reserve(response.n_elem);
    for (const double value : response) {
      decomposed_response.push_back(decompose_binary64(value));
    }
    const ExactDotPlan information_plan = make_exact_dot_plan(
      design_range_,
      binary_range(decomposed_response),
      rows_
    );
    std::vector<HighPrecision> information(columns_, HighPrecision(0));
    for (std::size_t column = 0; column < columns_; ++column) {
      information[column] = exact_dot(
        design_.data() + rows_ * column,
        decomposed_response.data(),
        rows_,
        information_plan
      );
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
  std::vector<Binary64> design_;
  BinaryRange design_range_;
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
