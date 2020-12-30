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

TEST(SphereTetOverlap, EdgeCases) {
	const std::vector<Sphere> spheres = {
		{vector_t::Zero(), 1.0},
		{{-0.01725, 0, 0}, 1.0},
		{vector_t::Zero(), 1.0}
	};

	const std::vector<Tetrahedron> tetrahedra = {
		{{{
			{0.0357829, 0, 1.01271},
			{0, 0, 1.01271},
			{0.0356948, 0.0386086, 0.962075},
			{0, 0, 0.962075}
		}}},
		{{{
			{0.9667906976744187, 0, 3.098296812907414e-16},
			{1.002654107311333, 0.0384643285369352, -2.82302142880589e-16},
			{1.002573643410853, 0, 4.131062417209885e-16},
			{1.002573643410853, 0, -0.05063534883720874}
		}}},
		{{{
			{0.28, -0.9599999999999999, -0.02102000000000028},
			{0.2400000000000001, -0.9599999999999999, 0.01898000000000015},
			{0.28, -0.9999999999999999, 0.01898000000000015},
			{0.28, -0.9599999999999999, 0.01898000000000015}
		}}}
	};

	ASSERT_EQ(spheres.size(), tetrahedra.size());

	const scalar_t epsilon =
		std::sqrt(std::numeric_limits<scalar_t>::epsilon());

	for(std::size_t idx = 0; idx < spheres.size(); ++idx) {
		std::array<Tetrahedron, 4> tets4;
		decompose(tetrahedra[idx], tets4);

		const scalar_t overlapCalcTet = overlap(spheres[idx], tetrahedra[idx]);

		scalar_t overlapCalcTets4 = 0;
		for(const auto& tet : tets4)
			overlapCalcTets4 += overlap(spheres[idx], tet);

		std::cout << "volume tet:  " << overlapCalcTet << std::endl;
		std::cout << "volume tets4:  " << overlapCalcTets4 << std::endl;

		ASSERT_NEAR(overlapCalcTet, overlapCalcTets4,
			epsilon * spheres[idx].volume);
	}
}
