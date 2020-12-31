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

#include "gtest/gtest.h"

#include "overlap.hpp"

#include "common.hpp"

TEST(DecomposeElements, Hexahedron) {
	Hexahedron hex = unitHexahedron();

	std::array<Tetrahedron, 4> subTets;
	std::array<Tetrahedron, 5> tets5;
	std::array<Tetrahedron, 6> tets6;
	std::array<Wedge, 2> wedges;

	decompose(hex, tets5);
	decompose(hex, tets6);
	decompose(hex, wedges);

	scalar_t tets5Volume = 0;
	for(const auto& tet : tets5)
		tets5Volume += tet.volume;

	scalar_t tets6Volume = 0;
	scalar_t tets24Volume = 0;
	for(const auto& tet : tets6) {
		decompose(tet, subTets);
		tets6Volume += tet.volume;

		for(const auto& subTet : subTets)
			tets24Volume += subTet.volume;
	}

	constexpr scalar_t delta(5e2 * std::numeric_limits<scalar_t>::epsilon());

	ASSERT_NEAR(hex.volume, tets5Volume, delta);
	ASSERT_NEAR(hex.volume, tets6Volume, delta);
	ASSERT_NEAR(hex.volume, tets24Volume, delta);
	ASSERT_NEAR(hex.volume, wedges[0].volume + wedges[1].volume, delta);
}
