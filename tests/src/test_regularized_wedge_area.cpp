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

// Test regularized_wedge_area using different values for the distance of the
// intersection point from the center of the sphere `z`.
TEST(RegularizedWedgeArea, Distance) {
  using namespace overlap::detail;

  // special cases should return precisely zero
  ASSERT_EQ(regularized_wedge_area(1.0, 1.0, 0.25 * pi), 0.0);
  ASSERT_EQ(regularized_wedge_area(1.0, -1.0, 0.25 * pi), 0.0);

  constexpr auto delta = Scalar{2e2} * std::numeric_limits<Scalar>::epsilon();

  ASSERT_NEAR(regularized_wedge_area(1.0, tinyEpsilon, 0.5 * pi), pi,
              5 * delta);

  ASSERT_NEAR(regularized_wedge_area(1.0, -tinyEpsilon, 0.5 * pi), pi,
              5 * delta);
}

// Test regularized_wedge_area using different values of the angle `alpha`.
TEST(regularized_wedge_area, Angle) {
  using namespace overlap::detail;

  // special cases should return constants values
  ASSERT_EQ(regularized_wedge_area(1.0, 0.0, 0.0), 0.0);
  ASSERT_EQ(regularized_wedge_area(1.0, 0.0, 0.5 * pi), pi);

  ASSERT_NEAR(regularized_wedge_area(1.0, 0.0, 0.75 * pi),
              2 * pi - regularized_wedge_area(1.0, 0.0, 0.25 * pi),
              std::numeric_limits<Scalar>::epsilon());
}
