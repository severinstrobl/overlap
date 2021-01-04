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

// Test regularizedWedge using different values of the distance `d`.
TEST(RegularizedWedge, Distance) {
	using namespace detail;

	// special case should return precisely zero
	ASSERT_EQ(regularizedWedge(1.0, 1.0, 0.25 * M_PI), 0.0);

	constexpr scalar_t epsilon = 5 * std::numeric_limits<scalar_t>::epsilon();

	ASSERT_NEAR(regularizedWedge(1.0, detail::tinyEpsilon, 0.25 * M_PI),
		M_PI / 6.0, epsilon);

	ASSERT_NEAR(regularizedWedge(1.0, detail::tinyEpsilon, 0.5 * M_PI),
		M_PI / 3.0, epsilon);
}

// Test regularizedWedge using different values of the angle `alpha`.
TEST(RegularizedWedge, Angle) {
	using namespace detail;

	// special case should return precisely zero
	ASSERT_EQ(regularizedWedge(1.0, 0.5, 0.0), 0.0);

	constexpr scalar_t epsilon = std::numeric_limits<scalar_t>::epsilon();

	ASSERT_NEAR(regularizedWedge(1.0, 0.5, 0.5 * M_PI), 5.0 * M_PI / 48.0,
		epsilon);

	// test using angle of alpha = pi/2
	constexpr scalar_t alpha = 0.5 * M_PI;
	constexpr scalar_t delta = std::numeric_limits<scalar_t>::epsilon();

	// introduce slight variations to `alpha` and `z`
	ASSERT_NEAR(regularizedWedge(1.0, 0.5, alpha,
			0.5 * std::cos(alpha + 0.5 * M_PI)),
		regularizedWedge(1.0, 0.5, alpha - delta,
			0.5 * std::cos(alpha + 0.5 * M_PI - delta)),
		5 * epsilon);

	ASSERT_NEAR(regularizedWedge(1.0, 0.5, alpha,
			0.5 * std::cos(alpha + 0.5 * M_PI)),
		regularizedWedge(1.0, 0.5, alpha + delta,
			0.5 * std::cos(alpha + 0.5 * M_PI + delta)),
		5 * epsilon);

	ASSERT_NEAR(regularizedWedge(1.0, 0.5, alpha,
			-0.5 * std::cos(alpha + 0.5 * M_PI)),
		regularizedWedge(1.0, 0.5, alpha - delta,
			-0.5 * std::cos(alpha + 0.5 * M_PI - delta)),
		5 * epsilon);

	ASSERT_NEAR(regularizedWedge(1.0, 0.5, alpha,
			-0.5 * std::cos(alpha + 0.5 * M_PI)),
		regularizedWedge(1.0, 0.5, alpha + delta,
			-0.5 * std::cos(alpha + 0.5 * M_PI + delta)),
		5 * epsilon);
}

// Test regularizedWedge for simple angles and base points very close to the
// center.
TEST(RegularizedWedge, NearCenter) {
	using namespace detail;

	constexpr scalar_t epsilon = 5 * std::numeric_limits<scalar_t>::epsilon();

	ASSERT_NEAR(regularizedWedge(1.0,
		std::numeric_limits<scalar_t>::epsilon(), 0.25 * M_PI), M_PI / 6.0,
		epsilon);

	ASSERT_NEAR(regularizedWedge(1.0,
		std::numeric_limits<scalar_t>::epsilon(), 0.5 * M_PI), M_PI / 3.0,
		epsilon);
}
