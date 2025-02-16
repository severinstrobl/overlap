// Copyright (C) 2025 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("DoublePrecision") {
  using namespace overlap;

  TEST_CASE("ConstantFloat") {
    // 32 bit floating point: 2^(24 - 24/2) + 1 = 2^12 + 1 = 4097
    CHECK_EQ(detail::DoublePrecision<float>::constant(), 4097);
  }

  TEST_CASE("ConstantDouble") {
    // 64 bit floating point: 2^(53 - int(53/2)) + 1 = 2^27 + 1 = 134217729
    CHECK_EQ(detail::DoublePrecision<double>::constant(), 134217729);
  }

  TEST_CASE_TEMPLATE("Splitting", T, float, double) {
    const auto value = static_cast<T>(detail::pi);
    const auto [high, low] = detail::DoublePrecision<T>::split(value);

    CHECK_GT(high, low);
    CHECK_EQ(high + low, value);

#if !defined(EIGEN_HAS_SINGLE_INSTRUCTION_MADD)
    auto [eigen_high, eigen_low] = std::pair{T{0}, T{0}};
    Eigen::internal::veltkamp_splitting(value, eigen_high, eigen_low);

    CHECK_EQ(high, eigen_high);
    CHECK_EQ(low, eigen_low);
#endif
  }

  TEST_CASE_TEMPLATE("TwoProduct", T, float, double) {
    const auto a = static_cast<T>(detail::pi);
    const auto b = std::numeric_limits<T>::epsilon() * a;

    const auto result = detail::DoublePrecision<T>::two_product(a, b);
    CHECK_EQ(result.high(), a * b);

    if constexpr (std::is_same_v<T, float>) {
      CHECK_EQ(static_cast<double>(result.high()) +
                   static_cast<double>(result.low()),
               static_cast<double>(a) * static_cast<double>(b));
    }

#if !defined(EIGEN_HAS_SINGLE_INSTRUCTION_MADD)
    auto [eigen_high, eigen_low] = std::pair{T{0}, T{0}};
    Eigen::internal::twoprod(a, b, eigen_high, eigen_low);

    CHECK_EQ(result.high(), eigen_high);
    CHECK_EQ(result.low(), eigen_low);
#endif
  }

  TEST_CASE_TEMPLATE("Add", T, float, double) {
    using DoublePrecision = detail::DoublePrecision<T>;

    const auto value = static_cast<T>(detail::pi);
    const auto delta = std::numeric_limits<T>::epsilon();

    CHECK_NE((value + delta) - value, delta);

    const auto result = (DoublePrecision{value} + DoublePrecision{delta}) +
                        DoublePrecision{-value};

    CHECK_EQ(result.value(), delta);

    if constexpr (std::is_same_v<T, float>) {
      CHECK_EQ(result.template as<double>(),
               (double{value} + double{delta}) - double{value});
    }
  }

  TEST_CASE_TEMPLATE("ConstExprAdd", T, float, double) {
    using DoublePrecision = detail::DoublePrecision<T>;

    constexpr auto value = static_cast<T>(EIGEN_PI);
    constexpr auto delta = std::numeric_limits<T>::epsilon();

    static_assert((value + delta) - value != delta);

    constexpr auto result = (DoublePrecision{value} + DoublePrecision{delta}) +
                            DoublePrecision{-value};

    static_assert(result.value() == delta);

    if constexpr (std::is_same_v<T, float>) {
      static_assert(result.template as<double>() ==
                    (double{value} + double{delta}) - double{value});
    }
  }

  TEST_CASE_TEMPLATE("Sub", T, float, double) {
    using DoublePrecision = detail::DoublePrecision<T>;

    const auto value = static_cast<T>(detail::pi);
    const auto delta = std::numeric_limits<T>::epsilon();

    CHECK_NE((value - delta) - value, -delta);

    const auto result = (DoublePrecision{value} - DoublePrecision{delta}) -
                        DoublePrecision{value};

    CHECK_EQ(result.value(), -delta);

    if constexpr (std::is_same_v<T, float>) {
      CHECK_EQ(result.template as<double>(),
               (double{value} - double{delta}) - double{value});
    }
  }

  TEST_CASE_TEMPLATE("AddSub", T, float, double) {
    using DoublePrecision = detail::DoublePrecision<T>;

    const auto value = static_cast<T>(detail::pi);
    const auto delta = std::numeric_limits<T>::epsilon();

    CHECK_NE((value + delta) - value, delta);

    const auto result = (DoublePrecision{value} + DoublePrecision{delta}) -
                        DoublePrecision{value};

    CHECK_EQ(result.value(), delta);
  }

  TEST_CASE_TEMPLATE("ConstExprAddSub", T, float, double) {
    using DoublePrecision = detail::DoublePrecision<T>;

    constexpr auto value = static_cast<T>(EIGEN_PI);
    constexpr auto delta = std::numeric_limits<T>::epsilon();

    static_assert((value + delta) - value != delta);

    constexpr auto result = (DoublePrecision{value} + DoublePrecision{delta}) -
                            DoublePrecision{value};

    static_assert(result.value() == delta);

    if constexpr (std::is_same_v<T, float>) {
      static_assert(result.template as<double>() ==
                    (double{value} + double{delta}) - double{value});
    }
  }

  TEST_CASE_TEMPLATE("Multiply", T, float, double) {
    using DoublePrecision = detail::DoublePrecision<T>;

    const auto a = static_cast<T>(detail::pi);
    const auto b = std::numeric_limits<T>::epsilon() * a;
    const auto c = std::sqrt(a);
    const auto d = DoublePrecision::two_product(a, b);
    const auto e = DoublePrecision::two_product(b, c);

    const auto result = d * e;
    CHECK_NE(result.value(), T{});

    if constexpr (std::is_same_v<T, float>) {
      const auto error =
          overlap::abs(result.template as<double>() -
                       (double{a} * double{b}) * (double{b} * double{c}));

      CHECK_LT(error, 8.0 * std::numeric_limits<double>::epsilon() *
                          result.template as<double>());
    }

#if !defined(EIGEN_HAS_SINGLE_INSTRUCTION_MADD)
    auto [eigen_high, eigen_low] = std::pair{T{0}, T{0}};
    Eigen::internal::twoprod(d.high(), d.low(), e.high(), e.low(), eigen_high,
                             eigen_low);

    CHECK_EQ(result.high(), eigen_high);
    CHECK_EQ(result.low(), eigen_low);
#endif
  }

  TEST_CASE_TEMPLATE("ConstExprMultiply", T, float, double) {
    using DoublePrecision = detail::DoublePrecision<T>;

    constexpr auto a = static_cast<T>(EIGEN_PI);
    constexpr auto b = std::numeric_limits<T>::epsilon() * a;
    constexpr auto c = a / 3.0;
    constexpr auto d = DoublePrecision::two_product(a, b);
    constexpr auto e = DoublePrecision::two_product(b, c);

    constexpr auto result = d * e;
    static_assert(result.value() != T{});

    if constexpr (std::is_same_v<T, float>) {
      constexpr auto error =
          overlap::abs(result.template as<double>() -
                       (double{a} * double{b}) * (double{b} * double{c}));

      static_assert(error < 4.0 * std::numeric_limits<double>::epsilon() *
                                result.template as<double>());
    }
  }
}
