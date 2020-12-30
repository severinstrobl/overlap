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

#ifndef OVERLAP_TEST_COMMON_HPP
#define OVERLAP_TEST_COMMON_HPP

#include <iostream>

#include "gtest/gtest.h"

#include "overlap.hpp"

inline Hexahedron unitHexahedron(scalar_t scaling = scalar_t(1)) {
	vector_t v0{-1, -1, -1};
	vector_t v1{ 1, -1, -1};
	vector_t v2{ 1,  1, -1};
	vector_t v3{-1,  1, -1};
	vector_t v4{-1, -1,  1};
	vector_t v5{ 1, -1,  1};
	vector_t v6{ 1,  1,  1};
	vector_t v7{-1,  1,  1};

	return Hexahedron{v0 * scaling, v1 * scaling, v2 * scaling, v3 * scaling,
		v4 * scaling, v5 * scaling, v6 * scaling, v7 * scaling};
}

inline void overlap(const Sphere& s, const Hexahedron& hex,
	const scalar_t epsilon, const scalar_t exactResult = scalar_t{-1}) {

	std::array<Tetrahedron, 4> subTets;
	std::array<Tetrahedron, 5> tets5;
	std::array<Tetrahedron, 6> tets6;
	std::array<Wedge, 2> wedges;

	decompose(hex, tets5);
	decompose(hex, tets6);
	decompose(hex, wedges);

	const scalar_t overlapCalcHex = overlap(s, hex);

	if(exactResult != scalar_t{-1})
        ASSERT_NEAR(overlapCalcHex, exactResult, epsilon);

	scalar_t overlapCalcTets5 = 0;
	for(const auto& tet : tets5)
		overlapCalcTets5 += overlap(s, tet);

	scalar_t overlapCalcTets6 = 0;
	scalar_t overlapCalcTets24 = 0;
	for(const auto& tet : tets6) {
		decompose(tet, subTets);

		for(const auto& subTet : subTets)
			overlapCalcTets24 += overlap(s, subTet);
			overlapCalcTets6 += overlap(s, tet);
	}

	scalar_t overlapCalcWedges = 0;
	for(const auto& wedge : wedges)
		overlapCalcWedges += overlap(s, wedge);

	std::cout << "volume hex:    " << overlapCalcHex << std::endl;
	std::cout << "volume wedges: " << overlapCalcWedges << std::endl;
	std::cout << "volume tets5:  " << overlapCalcTets5 << std::endl;
	std::cout << "volume tets6:  " << overlapCalcTets6 << std::endl;
	std::cout << "volume tets24: " << overlapCalcTets24 << std::endl;

	ASSERT_NEAR(overlapCalcHex, overlapCalcWedges, epsilon);
	ASSERT_NEAR(overlapCalcHex, overlapCalcTets5, epsilon);
	ASSERT_NEAR(overlapCalcHex, overlapCalcTets6, epsilon);
	ASSERT_NEAR(overlapCalcHex, overlapCalcTets24, epsilon);
	ASSERT_NEAR(overlapCalcTets5, overlapCalcTets6, epsilon);
	ASSERT_NEAR(overlapCalcTets5, overlapCalcTets24, epsilon);
	ASSERT_NEAR(overlapCalcTets6, overlapCalcTets24, epsilon);
}

inline void area(const Sphere& s, const Hexahedron& hex,
	const scalar_t epsilon) {

	std::array<Tetrahedron, 4> subTets;
	std::array<Tetrahedron, 5> tets5;
	std::array<Tetrahedron, 6> tets6;
	std::array<Wedge, 2> wedges;

	decompose(hex, tets5);
	decompose(hex, tets6);
	decompose(hex, wedges);

	const scalar_t areaCalcHex = overlapArea(s, hex)[0];

	scalar_t areaCalcWedges = 0;
	scalar_t areaCalcTets5 = 0;
	scalar_t areaCalcTets6 = 0;
	scalar_t areaCalcTets24 = 0;

	for(const auto& tet : tets5)\
		areaCalcTets5 += overlapArea(s, tet)[0];

	for(const auto& tet : tets6) {
		decompose(tet, subTets);

		for(const auto& subTet : subTets)
			areaCalcTets24 += overlapArea(s, subTet)[0];
			areaCalcTets6 += overlapArea(s, tet)[0];
	}

	for(const auto& wedge : wedges)
		areaCalcWedges += overlapArea(s, wedge)[0];

	std::cout << "sphere center: [" << s.center.transpose() << "], radius: " <<
		s.radius << std::endl;

	std::cout << "sphere surface hex:    " << areaCalcHex << std::endl;
	std::cout << "sphere surface wedges: " << areaCalcWedges << std::endl;
	std::cout << "sphere surface tets5:  " << areaCalcTets5 << std::endl;
	std::cout << "sphere surface tets6:  " << areaCalcTets6 << std::endl;
	std::cout << "sphere surface tets24: " << areaCalcTets24 << std::endl;

	ASSERT_NEAR(areaCalcHex, areaCalcWedges, epsilon);
	ASSERT_NEAR(areaCalcHex, areaCalcTets5, epsilon);
	ASSERT_NEAR(areaCalcHex, areaCalcTets6, epsilon);
	ASSERT_NEAR(areaCalcHex, areaCalcTets24, epsilon);
	ASSERT_NEAR(areaCalcTets5, areaCalcTets6, epsilon);
	ASSERT_NEAR(areaCalcTets5, areaCalcTets24, epsilon);
	ASSERT_NEAR(areaCalcTets6, areaCalcTets24, epsilon);
}

#endif // OVERLAP_TEST_COMMON_HPP
