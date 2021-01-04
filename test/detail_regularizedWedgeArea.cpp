/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2015-2020 Severin Strobl <severin.strobl@fau.de>
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

#define _USE_MATH_DEFINES

#include "gtest/gtest.h"

#include "overlap.hpp"

// Test regularizedWedgeArea using different values for the distance of the
// intersection point from the center of the sphere `z`.
TEST(RegularizedWedgeArea, Distance) {
	using namespace detail;

	// special cases should return precisely zero
	ASSERT_EQ(regularizedWedgeArea(1.0, 1.0, 0.25 * M_PI), 0.0);
	ASSERT_EQ(regularizedWedgeArea(1.0, -1.0, 0.25 * M_PI), 0.0);

	constexpr scalar_t delta(2e2 * std::numeric_limits<scalar_t>::epsilon());

	ASSERT_NEAR(regularizedWedgeArea(1.0, detail::tinyEpsilon,
		0.5 * M_PI), M_PI, 5 * delta);

	ASSERT_NEAR(regularizedWedgeArea(1.0, -detail::tinyEpsilon,
		0.5 * M_PI), M_PI, 5 * delta);
}

// Test regularizedWedgeArea using different values of the angle `alpha`.
TEST(regularizedWedgeArea, Angle) {
	using namespace detail;

	// special cases should return constants values
	ASSERT_EQ(regularizedWedgeArea(1.0, 0.0, 0.0), 0.0);
	ASSERT_EQ(regularizedWedgeArea(1.0, 0.0, 0.5 * M_PI), M_PI);

	ASSERT_NEAR(regularizedWedgeArea(1.0, 0.0, 0.75 * M_PI),
		2 * M_PI - regularizedWedgeArea(1.0, 0.0, 0.25 * M_PI),
		std::numeric_limits<scalar_t>::epsilon());
}
