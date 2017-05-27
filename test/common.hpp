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

#ifndef COMMON_HPP
#define COMMON_HPP

#include <iostream>

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

#define overlapWorker(s, hex, delta, exactResult) {\
	std::array<Tetrahedron, 4> subTets;\
	std::array<Tetrahedron, 5> tets5;\
	std::array<Tetrahedron, 6> tets6;\
	std::array<Wedge, 2> wedges;\
\
	decompose(hex, tets5);\
	decompose(hex, tets6);\
	decompose(hex, wedges);\
\
	scalar_t overlapCalcHex = overlap(s, hex);\
\
	if(exactResult != scalar_t(-1)) {\
        if(std::abs(exactResult - delta) <= delta)\
            BOOST_CHECK_SMALL(overlapCalcHex, delta);\
        else\
            BOOST_CHECK_CLOSE(overlapCalcHex, exactResult, delta);\
    }\
\
	scalar_t overlapCalcWedges = 0;\
	scalar_t overlapCalcTets5 = 0;\
	scalar_t overlapCalcTets6 = 0;\
	scalar_t overlapCalcTets24 = 0;\
\
	for(const auto& tet : tets5)\
		overlapCalcTets5 += overlap(s, tet);\
\
	for(const auto& tet : tets6) {\
		decompose(tet, subTets);\
\
		for(const auto& subTet : subTets)\
			overlapCalcTets24 += overlap(s, subTet);\
			overlapCalcTets6 += overlap(s, tet);\
	}\
\
	for(const auto& wedge : wedges)\
		overlapCalcWedges += overlap(s, wedge);\
\
	std::cout << "volume hex:    " << overlapCalcHex << std::endl;\
	std::cout << "volume wedges: " << overlapCalcWedges << std::endl;\
	std::cout << "volume tets5:  " << overlapCalcTets5 << std::endl;\
	std::cout << "volume tets6:  " << overlapCalcTets6 << std::endl;\
	std::cout << "volume tets24: " << overlapCalcTets24 << std::endl;\
\
	BOOST_CHECK_CLOSE(overlapCalcHex, overlapCalcWedges, delta);\
	BOOST_CHECK_CLOSE(overlapCalcHex, overlapCalcTets5, delta);\
	BOOST_CHECK_CLOSE(overlapCalcHex, overlapCalcTets6, delta);\
	BOOST_CHECK_CLOSE(overlapCalcHex, overlapCalcTets24, delta);\
	BOOST_CHECK_CLOSE(overlapCalcTets5, overlapCalcTets6, delta);\
	BOOST_CHECK_CLOSE(overlapCalcTets5, overlapCalcTets24, delta);\
	BOOST_CHECK_CLOSE(overlapCalcTets6, overlapCalcTets24, delta);\
}\

#endif // COMMON_HPP
