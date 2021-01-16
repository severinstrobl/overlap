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

#include "gtest/gtest.h"

#include "overlap.hpp"

TEST(Wedge, Constructors) {
	Wedge wedge0{};

	ASSERT_EQ(wedge0.volume, 0.0);

	Wedge wedge1{
		vector_t{-1, -1, -1}, vector_t{1, -1, -1}, vector_t{1,  1, -1},
		vector_t{-1, -1,  1}, vector_t{1, -1,  1}, vector_t{1,  1,  1}
	};

	Wedge wedge2{std::array<vector_t, 6>{{
		{-1, -1, -1}, {1, -1, -1}, {1,  1, -1},
		{-1, -1,  1}, {1, -1,  1}, {1,  1,  1}
	}}};

	constexpr scalar_t eps = std::numeric_limits<scalar_t>::epsilon();

	ASSERT_NEAR(wedge1.volume, wedge2.volume, eps);
}

TEST(Hexahedron, Constructors) {
	Hexahedron hex1{
		vector_t{-1, -1, -1}, vector_t{1, -1, -1},
		vector_t{ 1,  1, -1}, vector_t{-1,  1, -1},

		vector_t{-1, -1,  1}, vector_t{1, -1,  1},
		vector_t{ 1,  1,  1}, vector_t{-1,  1,  1}
	};

	Hexahedron hex2{std::array<vector_t, 8>{{
		{-1, -1, -1}, {1, -1, -1}, {1,  1, -1}, {-1,  1, -1},
		{-1, -1,  1}, {1, -1,  1}, {1,  1,  1}, {-1,  1,  1}
	}}};

	constexpr scalar_t eps = std::numeric_limits<scalar_t>::epsilon();

	ASSERT_NEAR(hex1.volume, hex2.volume, eps);
}
