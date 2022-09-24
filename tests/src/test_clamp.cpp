/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2021-2022 Severin Strobl <git@severin-strobl.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

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

  // Test clamping of numbers at lower/upper limit.
  TEST_CASE("Limits") {
    // clang-format off
    CHECK_EQ(overlap::detail::clamp(-1.0, -1.0, 1.0, 0.0), -1.0);
    CHECK_EQ(overlap::detail::clamp( 1.0, -1.0, 1.0, 0.0),  1.0);
    // clang-format on
  }

  // Test error handling of clamping of numbers using invalid inputs in debug
  // mode.
  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
  TEST_CASE("InvalidArguments") {
    // clang-format off
    REQUIRE_THROWS_AS(overlap::detail::clamp(0.0,  1.0, -1.0), AssertionError);
    REQUIRE_THROWS_AS(overlap::detail::clamp(0.0, -1.0,  1.0, -1.0), AssertionError);
    // clang-format on
  }
}