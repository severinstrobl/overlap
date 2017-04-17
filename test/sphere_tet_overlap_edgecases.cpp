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

#define BOOST_TEST_MODULE sphere_tet_overlap_edgecases
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "overlap.hpp"

#include "common.hpp"

BOOST_AUTO_TEST_CASE(sphere_tet_overlap_edgecase_0) {
	Sphere s(vector_t::Zero(), scalar_t(1));

	vector_t v0{0.0357829, 0, 1.01271};
	vector_t v1{0, 0, 1.01271};
	vector_t v2{0.0356948, 0.0386086, 0.962075};
	vector_t v3{0, 0, 0.962075};

	Tetrahedron tet = {v0, v1, v2, v3};

	scalar_t delta = std::sqrt(std::numeric_limits<scalar_t>::epsilon()) *
		s.volume;

	std::array<Tetrahedron, 4> tets4;
	decompose(tet, tets4);

	scalar_t overlapCalcTet = overlap(s, tet);
	scalar_t overlapCalcTets4 = 0;

	for(const auto& tet : tets4)
		overlapCalcTets4 += overlap(s, tet);

  	std::cout << "volume tet:  " << overlapCalcTet << std::endl;
	std::cout << "volume tets4:  " << overlapCalcTets4 << std::endl;

	BOOST_CHECK_CLOSE(overlapCalcTet, overlapCalcTets4, delta);
}

BOOST_AUTO_TEST_CASE(sphere_tet_overlap_edgecase_1) {
	Sphere s(vector_t{-0.01725, 0, 0}, scalar_t(1));

	vector_t v0{0.9667906976744187, 0, 3.098296812907414e-16};
	vector_t v1{1.002654107311333, 0.0384643285369352, -2.82302142880589e-16};
	vector_t v2{1.002573643410853, 0, 4.131062417209885e-16};
	vector_t v3{1.002573643410853, 0, -0.05063534883720874};

	Tetrahedron tet = {v0, v1, v2, v3};

	scalar_t delta = std::sqrt(std::numeric_limits<scalar_t>::epsilon()) *
		s.volume;

	std::array<Tetrahedron, 4> tets4;
	decompose(tet, tets4);

	scalar_t overlapCalcTet = overlap(s, tet);
	scalar_t overlapCalcTets4 = 0;

	for(const auto& tet : tets4)
		overlapCalcTets4 += overlap(s, tet);

  	std::cout << "volume tet:  " << overlapCalcTet << std::endl;
	std::cout << "volume tets4:  " << overlapCalcTets4 << std::endl;

	BOOST_CHECK_CLOSE(overlapCalcTet, overlapCalcTets4, delta);
}

BOOST_AUTO_TEST_CASE(sphere_tet_overlap_edgecase_2) {
	Sphere s(vector_t::Zero(), scalar_t(1));

	vector_t v0{0.28, -0.9599999999999999, -0.02102000000000028};
	vector_t v1{0.2400000000000001, -0.9599999999999999, 0.01898000000000015};
	vector_t v2{0.28, -0.9999999999999999, 0.01898000000000015};
	vector_t v3{0.28, -0.9599999999999999, 0.01898000000000015};

	Tetrahedron tet = {v0, v1, v2, v3};

	scalar_t delta = std::sqrt(std::numeric_limits<scalar_t>::epsilon()) *
		s.volume;

	std::array<Tetrahedron, 4> tets4;
	decompose(tet, tets4);

	scalar_t overlapCalcTet = overlap(s, tet);
	scalar_t overlapCalcTets4 = 0;

	for(const auto& tet : tets4)
		overlapCalcTets4 += overlap(s, tet);

  	std::cout << "volume tet:  " << overlapCalcTet << std::endl;
	std::cout << "volume tets4:  " << overlapCalcTets4 << std::endl;

	BOOST_CHECK_CLOSE(overlapCalcTet, overlapCalcTets4, delta);
}
