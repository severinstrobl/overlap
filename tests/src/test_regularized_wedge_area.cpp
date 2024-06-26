// Copyright (C) 2015-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("RegularizedWedgeArea") {
  using namespace overlap::detail;

  // Test regularized_wedge_area using different values for the distance of the
  // intersection point from the center of the sphere `z`.
  TEST_CASE("Distance") {
    // special cases should return precisely zero
    CHECK(regularized_wedge_area(1.0, 1.0, 0.25 * pi) == Scalar{0});
    CHECK(regularized_wedge_area(1.0, -1.0, 0.25 * pi) == Scalar{0});

    constexpr auto epsilon =
        Scalar{1e3} * std::numeric_limits<Scalar>::epsilon();

    CHECK(regularized_wedge_area(1.0, tiny_epsilon, 0.5 * pi) ==
          Approx(pi).epsilon(epsilon));

    CHECK(regularized_wedge_area(1.0, -tiny_epsilon, 0.5 * pi) ==
          Approx(pi).epsilon(epsilon));
  }

  // Test regularized_wedge_area using different values of the angle `alpha`.
  TEST_CASE("Angle") {
    // special cases should return constants values
    CHECK(regularized_wedge_area(1.0, 0.0, 0.0) == Scalar{0});
    CHECK(regularized_wedge_area(1.0, 0.0, 0.5 * pi) == pi);

    CHECK(regularized_wedge_area(1.0, 0.0, 0.75 * pi) ==
          Approx(2 * pi - regularized_wedge_area(1.0, 0.0, 0.25 * pi))
              .epsilon(std::numeric_limits<Scalar>::epsilon()));
  }

#ifndef NDEBUG
  // Test the clamping of input arguments slightly out of range.
  TEST_CASE("TestClamping") {
    CHECK(regularized_wedge_area(
              1.0, 0.5, -std::numeric_limits<Scalar>::epsilon()) == Scalar{0});

    CHECK(regularized_wedge_area(
              1.0, 0.0, pi + 1.5 * std::numeric_limits<Scalar>::epsilon()) ==
          Approx(2 * pi).epsilon(tiny_epsilon));
  }
#endif  // NDEBUG

  // Test error handling using invalid inputs in debug mode.
  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
  TEST_CASE("InvalidArguments") {
    REQUIRE_THROWS_AS(regularized_wedge_area(0.0, 1.0, 0.25 * pi),
                      AssertionError);

    REQUIRE_THROWS_AS(regularized_wedge_area(1, -2.0, 0.25 * pi),
                      AssertionError);
    REQUIRE_THROWS_AS(regularized_wedge_area(1, +2.0, 0.25 * pi),
                      AssertionError);

    REQUIRE_THROWS_AS(regularized_wedge_area(1, 0.5, -0.25 * pi),
                      AssertionError);
    REQUIRE_THROWS_AS(regularized_wedge_area(1, 0.5, 2 * pi), AssertionError);
  }
}
