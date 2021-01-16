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

TEST(Polygon, IsPlanar) {
	Triangle tri{
		vector_t{0.0, 0.0, 0.0},
		vector_t{1.0, 0.0, 0.0},
		vector_t{1.0, 1.0, 0.0}
	};

	ASSERT_TRUE(tri.isPlanar());

	Quadrilateral quad{
		vector_t{0.0, 0.0, 0.0},
		vector_t{1.0, 0.0, 0.0},
		vector_t{1.0, 1.0, 0.0},
		vector_t{0.0, 1.0, 0.01}
	};

	ASSERT_FALSE(quad.isPlanar());
}
