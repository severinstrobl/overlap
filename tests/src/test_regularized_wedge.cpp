/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2015-2022 Severin Strobl <severin.strobl@fau.de>
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

#include "gtest/gtest.h"

#include "overlap/overlap.hpp"

// Test regularized_wedge using different values of the distance `d`.
TEST(RegularizedWedge, Distance) {
  using namespace overlap::detail;

  // special case should return precisely zero
  ASSERT_EQ(regularized_wedge(1.0, 1.0, 0.25 * pi), 0.0);

  constexpr auto epsilon = 5 * std::numeric_limits<Scalar>::epsilon();

  ASSERT_NEAR(regularized_wedge(1.0, tinyEpsilon, 0.25 * pi), pi / 6.0, epsilon);
  ASSERT_NEAR(regularized_wedge(1.0, tinyEpsilon, 0.5 * pi), pi / 3.0, epsilon);
}

// Test regularized_wedge using different values of the angle `alpha`.
TEST(RegularizedWedge, Angle) {
  using namespace overlap::detail;

  // special case should return precisely zero
  ASSERT_EQ(regularized_wedge(1.0, 0.5, 0.0), 0.0);

  constexpr auto epsilon = std::numeric_limits<Scalar>::epsilon();

  ASSERT_NEAR(regularized_wedge(1.0, 0.5, 0.5 * pi), 5.0 * pi / 48.0, epsilon);

  // test using angle of alpha = pi/2
  const auto alpha = (Scalar{1} / Scalar{2}) * pi;
  constexpr auto delta = std::numeric_limits<Scalar>::epsilon();

  // introduce slight variations to `alpha` and `z`
  ASSERT_NEAR(
      regularized_wedge(1.0, 0.5, alpha, 0.5 * std::cos(alpha + 0.5 * pi)),
      regularized_wedge(1.0, 0.5, alpha - delta,
                       0.5 * std::cos(alpha + 0.5 * pi - delta)),
      5 * epsilon);

  ASSERT_NEAR(
      regularized_wedge(1.0, 0.5, alpha, 0.5 * std::cos(alpha + 0.5 * pi)),
      regularized_wedge(1.0, 0.5, alpha + delta,
                       0.5 * std::cos(alpha + 0.5 * pi + delta)),
      5 * epsilon);

  ASSERT_NEAR(
      regularized_wedge(1.0, 0.5, alpha, -0.5 * std::cos(alpha + 0.5 * pi)),
      regularized_wedge(1.0, 0.5, alpha - delta,
                       -0.5 * std::cos(alpha + 0.5 * pi - delta)),
      5 * epsilon);

  ASSERT_NEAR(
      regularized_wedge(1.0, 0.5, alpha, -0.5 * std::cos(alpha + 0.5 * pi)),
      regularized_wedge(1.0, 0.5, alpha + delta,
                       -0.5 * std::cos(alpha + 0.5 * pi + delta)),
      5 * epsilon);
}

// Test regularized_wedge for simple angles and base points very close to the
// center.
TEST(RegularizedWedge, NearCenter) {
  using namespace overlap::detail;

  constexpr auto epsilon = 5 * std::numeric_limits<Scalar>::epsilon();

  ASSERT_NEAR(
      regularized_wedge(1.0, std::numeric_limits<Scalar>::epsilon(), 0.25 * pi),
      pi / 6.0, epsilon);

  ASSERT_NEAR(
      regularized_wedge(1.0, std::numeric_limits<Scalar>::epsilon(), 0.5 * pi),
      pi / 3.0, epsilon);
}
