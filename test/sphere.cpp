/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2021 Severin Strobl <git@severin-strobl.de>
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

TEST(Sphere, Volume) {
	Sphere s{vector_t::Zero(), 1.0};


	constexpr scalar_t eps = std::numeric_limits<scalar_t>::epsilon();

	ASSERT_NEAR(s.volume, 4.0 / 3.0 * M_PI, eps);

	ASSERT_EQ(s.capVolume(-1 * s.radius), 0.0);
	ASSERT_EQ(s.capVolume(0), 0.0);
	ASSERT_NEAR(s.capVolume(0.5 * s.radius), (M_PI * 0.25 / 3.0) * 2.5, eps);
	ASSERT_NEAR(s.capVolume(s.radius), 0.5 * s.volume, eps);
	ASSERT_EQ(s.capVolume(2 * s.radius), s.volume);
	ASSERT_EQ(s.capVolume(3 * s.radius), s.volume);
}

TEST(Sphere, SurfaceArea) {
	Sphere s{vector_t::Zero(), 1.0};

	constexpr scalar_t eps = std::numeric_limits<scalar_t>::epsilon();

	ASSERT_NEAR(s.surfaceArea(), 4.0 * M_PI, eps);

	ASSERT_EQ(s.capSurfaceArea(-1 * s.radius), 0.0);
	ASSERT_EQ(s.capSurfaceArea(0), 0.0);
	ASSERT_EQ(s.capSurfaceArea(s.radius), 0.5 * s.surfaceArea());
	ASSERT_EQ(s.capSurfaceArea(2 * s.radius), s.surfaceArea());
	ASSERT_EQ(s.capSurfaceArea(3 * s.radius), s.surfaceArea());
}

TEST(Sphere, DiskArea) {
	Sphere s{vector_t::Zero(), 1.0};

	constexpr scalar_t eps = std::numeric_limits<scalar_t>::epsilon();

	ASSERT_EQ(s.diskArea(-1 * s.radius), 0.0);
	ASSERT_EQ(s.diskArea(0), 0.0);
	ASSERT_NEAR(s.diskArea(s.radius), M_PI, eps);
	ASSERT_EQ(s.diskArea(2 * s.radius), 0.0);
	ASSERT_EQ(s.diskArea(3 * s.radius), 0.0);
}

TEST(Sphere, Contains) {
	Sphere s{vector_t::Zero(), 2.0};

	ASSERT_TRUE(contains(s, vector_t{1, 1, 1}));
	ASSERT_FALSE(contains(s, vector_t{2, 2, 2}));
}
