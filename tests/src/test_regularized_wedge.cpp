// Copyright (C) 2015-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("RegularizedWedge") {
  using namespace overlap::detail;

  // Test regularized_wedge using different values of the distance `d`.
  TEST_CASE("VaryingDistance") {
    // special case should return precisely zero
    CHECK(regularized_wedge(1.0, 1.0, 0.25 * pi) == Scalar{0});

    constexpr auto epsilon = 5 * std::numeric_limits<Scalar>::epsilon();

    CHECK(regularized_wedge(1.0, tiny_epsilon, 0.25 * pi) ==
          Approx(pi / 6.0).epsilon(epsilon));
    CHECK(regularized_wedge(1.0, tiny_epsilon, 0.5 * pi) ==
          Approx(pi / 3.0).epsilon(epsilon));
  }

  constexpr auto epsilon = Scalar{5} * std::numeric_limits<Scalar>::epsilon();

  // Test regularized_wedge using different values of the angle `alpha`.
  TEST_CASE("VaryingAngle") {
    // special case should return precisely zero
    CHECK(regularized_wedge(1.0, 0.5, 0.0) == Scalar{0});

    CHECK(regularized_wedge(1.0, 0.5, 0.5 * pi) ==
          Approx(5.0 * pi / 48.0).epsilon(epsilon));
  }

  TEST_CASE("Stability") {
    // test using angle of alpha = pi/2
    const auto alpha = (Scalar{1} / Scalar{2}) * pi;
    constexpr auto delta = std::numeric_limits<Scalar>::epsilon();

    // introduce slight variations to `alpha` and `z`
    CHECK(
        regularized_wedge(1.0, 0.5, alpha, 0.5 * std::cos(alpha + 0.5 * pi)) ==
        Approx(regularized_wedge(1.0, 0.5, alpha - delta,
                                 0.5 * std::cos(alpha + 0.5 * pi - delta)))
            .epsilon(epsilon));

    CHECK(
        regularized_wedge(1.0, 0.5, alpha, 0.5 * std::cos(alpha + 0.5 * pi)) ==
        Approx(regularized_wedge(1.0, 0.5, alpha + delta,
                                 0.5 * std::cos(alpha + 0.5 * pi + delta)))
            .epsilon(epsilon));

    CHECK(
        regularized_wedge(1.0, 0.5, alpha, -0.5 * std::cos(alpha + 0.5 * pi)) ==
        Approx(regularized_wedge(1.0, 0.5, alpha - delta,
                                 -0.5 * std::cos(alpha + 0.5 * pi - delta)))
            .epsilon(epsilon));

    CHECK(
        regularized_wedge(1.0, 0.5, alpha, -0.5 * std::cos(alpha + 0.5 * pi)) ==
        Approx(regularized_wedge(1.0, 0.5, alpha + delta,
                                 -0.5 * std::cos(alpha + 0.5 * pi + delta)))
            .epsilon(epsilon));
  }

  // Test regularized_wedge for simple angles and base points very close to the
  // center.
  TEST_CASE("NearCenter") {
    CHECK(regularized_wedge(1.0, std::numeric_limits<Scalar>::epsilon(),
                            0.25 * pi) == Approx(pi / 6.0).epsilon(epsilon));

    CHECK(regularized_wedge(1.0, std::numeric_limits<Scalar>::epsilon(),
                            0.5 * pi) == Approx(pi / 3.0).epsilon(epsilon));
  }

#ifndef NDEBUG
  // Test the clamping of input arguments slightly out of range.
  TEST_CASE("TestClamping") {
    CHECK(regularized_wedge(
              1.0, 0.5, -std::numeric_limits<Scalar>::epsilon()) == Scalar{0});

    CHECK(
        regularized_wedge(1.0, tiny_epsilon,
                          0.5 * pi + std::numeric_limits<Scalar>::epsilon()) ==
        Approx(pi / 3.0).epsilon(epsilon));
  }
#endif  // NDEBUG

  // Test error handling using invalid inputs in debug mode.
  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
  TEST_CASE("InvalidArguments") {
    REQUIRE_THROWS_AS(regularized_wedge(0.0, 1.0, 0.25 * pi), AssertionError);
    REQUIRE_THROWS_AS(regularized_wedge(1, -1.0, 0.25 * pi), AssertionError);
    REQUIRE_THROWS_AS(regularized_wedge(1, 2.0, 0.25 * pi), AssertionError);
    REQUIRE_THROWS_AS(regularized_wedge(1, 0.5, -0.25 * pi), AssertionError);
    REQUIRE_THROWS_AS(regularized_wedge(1, 0.5, pi), AssertionError);
  }
}
