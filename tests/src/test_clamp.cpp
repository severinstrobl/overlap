// Copyright (C) 2021-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("clamp") {
  // Test clamping of numbers without tolerance.
  TEST_CASE("WithoutTolerance") {
    // clang-format off
    CHECK_EQ(overlap::detail::clamp(-0.5, -1.0, 1.0), -0.5);
    CHECK_EQ(overlap::detail::clamp( 0.5, -1.0, 1.0),  0.5);
    // clang-format on
  }

  // Test clamping of numbers with tolerance.
  TEST_CASE("WithTolerance") {
    // clang-format off
    CHECK_EQ(overlap::detail::clamp(-1.1, -1.0, 1.0, 0.25), -1.0);
    CHECK_EQ(overlap::detail::clamp( 1.1, -1.0, 1.0, 0.25),  1.0);
    // clang-format on
  }

  // Test clamping of numbers while exceeding the tolerance.
  TEST_CASE("WithToleranceExceeded") {
    // clang-format off
    CHECK_EQ(overlap::detail::clamp(-1.1, -1.0, 1.0, 0.01), -1.1);
    CHECK_EQ(overlap::detail::clamp( 1.1, -1.0, 1.0, 0.01),  1.1);
    // clang-format on
  }

  // Test clamping of numbers at lower/upper limit.
  TEST_CASE("Limits") {
    // clang-format off
    CHECK_EQ(overlap::detail::clamp(-1.0, -1.0, 1.0, 0.0), -1.0);
    CHECK_EQ(overlap::detail::clamp( 1.0, -1.0, 1.0, 0.0),  1.0);
    // clang-format on
  }

  // Test error handling of clamping of numbers using invalid inputs in debug
  // mode.
  TEST_CASE("InvalidArguments") {
    // clang-format off
    REQUIRE_THROWS_AS(overlap::detail::clamp(0.0,  1.0, -1.0), AssertionError);
    REQUIRE_THROWS_AS(overlap::detail::clamp(0.0, -1.0,  1.0, -1.0), AssertionError);
    REQUIRE_THROWS_AS(overlap::detail::clamp(0.0,  1.0, -1.0, -1.0), AssertionError);
    // clang-format on
  }
}
