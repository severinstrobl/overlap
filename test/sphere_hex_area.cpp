/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2015-2017 Severin Strobl <severin.strobl@fau.de>
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

#define BOOST_TEST_MODULE sphere_hex_area
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "overlap.hpp"

#include "common.hpp"

// Sphere intersects one face.
BOOST_AUTO_TEST_CASE(sphere_hex_area_face) {
	Hexahedron hex = unitHexahedron();
	Sphere s({0, 0, 1}, 0.75);

	auto result = overlapArea(s, hex);

	std::array<scalar_t, 8> resultExact;
	resultExact.fill(scalar_t(0));
	resultExact[0] = scalar_t(0.5) * s.surfaceArea();
	resultExact[6] = s.diskArea(s.radius);
	resultExact[7] = resultExact[6];

	scalar_t delta(1e2 * std::numeric_limits<scalar_t>::epsilon());
	for(size_t i = 0; i < resultExact.size(); ++i) {
		std::cout << "result[" << i << "]: " << result[i] << std::endl;
		BOOST_CHECK_CLOSE(result[i], resultExact[i], delta);
	}
}

// Sphere intersects one edge (and thus 1 edge and 2 faces).
BOOST_AUTO_TEST_CASE(sphere_hex_area_edge) {
	Hexahedron hex = unitHexahedron();
	Sphere s({1, 1, 0}, 0.75);

	auto result = overlapArea(s, hex);

	std::array<scalar_t, 8> resultExact;
	resultExact.fill(scalar_t(0));
	resultExact[0] = scalar_t(0.25) * s.surfaceArea();
	resultExact[3] = 0.5 * s.diskArea(s.radius);
	resultExact[4] = resultExact[3];
	resultExact[7] = 2 * resultExact[3];

	scalar_t delta(1e2 * std::numeric_limits<scalar_t>::epsilon());
	for(size_t i = 0; i < resultExact.size(); ++i) {
		std::cout << "result[" << i << "]: " << result[i] << std::endl;
		BOOST_CHECK_CLOSE(result[i], resultExact[i], delta);
	}
}

// Sphere intersects one vertex (and thus 3 edge and 3 faces).
BOOST_AUTO_TEST_CASE(sphere_hex_area_vertex) {
	Hexahedron hex = unitHexahedron();
	Sphere s({1, 1, 1}, 0.75);

	auto result = overlapArea(s, hex);

	std::array<scalar_t, 8> resultExact;
	resultExact.fill(scalar_t(0));
	resultExact[0] = scalar_t(0.125) * s.surfaceArea();
	resultExact[3] = 0.25 * s.diskArea(s.radius);
	resultExact[4] = resultExact[3];
	resultExact[6] = resultExact[3];
	resultExact[7] = 3 * resultExact[3];

	scalar_t delta(1e3 * std::numeric_limits<scalar_t>::epsilon());
	for(size_t i = 0; i < resultExact.size(); ++i) {
		std::cout << "result[" << i << "]: " << result[i] << std::endl;
		BOOST_CHECK_CLOSE(result[i], resultExact[i], delta);
	}
}
