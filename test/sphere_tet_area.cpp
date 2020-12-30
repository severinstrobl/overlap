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

#include <cmath>

#include "overlap.hpp"

#include "common.hpp"

// Sphere inside of element.
TEST(SphereTetAreaTest, SphereInTet) {
	vector_t v0{-std::sqrt(3) / 6.0, -1.0 / 2.0, 0};
	vector_t v1{std::sqrt(3) / 3.0, 0, 0};
	vector_t v2{-std::sqrt(3) / 6.0, +1.0 / 2.0, 0};
	vector_t v3{0, 0, std::sqrt(6) / 3.0};

	Tetrahedron tet{v0, v1, v2, v3};
	Sphere s({0, 0, 0.25}, 0.125);

	auto result = overlapArea(s, tet);

	std::array<scalar_t, 6> resultExact;
	resultExact.fill(scalar_t(0));
	resultExact[0] = s.surfaceArea();

	for(size_t i = 0; i < resultExact.size(); ++i)
		ASSERT_EQ(result[i], resultExact[i]);
}

// Element contained in sphere.
TEST(SphereTetAreaTest, TetInSphere) {
	vector_t v0{-std::sqrt(3) / 6.0, -1.0 / 2.0, 0};
	vector_t v1{std::sqrt(3) / 3.0, 0, 0};
	vector_t v2{-std::sqrt(3) / 6.0, +1.0 / 2.0, 0};
	vector_t v3{0, 0, std::sqrt(6) / 3.0};

	Tetrahedron tet{v0, v1, v2, v3};
	Sphere s({0, 0, std::sqrt(6) / 6.0}, 2);

	auto result = overlapArea(s, tet);

	std::array<scalar_t, 6> resultExact;
	resultExact[0] = scalar_t(0);
	resultExact[1] = tet.faces[0].area;
	resultExact[2] = tet.faces[1].area;
	resultExact[3] = tet.faces[2].area;
	resultExact[4] = tet.faces[3].area;
	resultExact.back() = std::accumulate(resultExact.begin() + 1,
		resultExact.end() - 1, scalar_t(0));

	for(size_t i = 0; i < resultExact.size(); ++i)
		ASSERT_NEAR(result[i], resultExact[i],
			std::numeric_limits<scalar_t>::epsilon());
}

// Sphere intersects one face.
TEST(SphereTetAreaTest, Face) {
	vector_t v0{-std::sqrt(3) / 6.0, -1.0 / 2.0, 0};
	vector_t v1{std::sqrt(3) / 3.0, 0, 0};
	vector_t v2{-std::sqrt(3) / 6.0, +1.0 / 2.0, 0};
	vector_t v3{0, 0, std::sqrt(6) / 3.0};

	Tetrahedron tet{v0, v1, v2, v3};
	Sphere s({0, 0, 0}, 0.25);

	auto result = overlapArea(s, tet);

	std::array<scalar_t, 6> resultExact;
	resultExact.fill(scalar_t(0));
	resultExact[0] = scalar_t(0.5) * s.surfaceArea();
	resultExact[1] = s.diskArea(s.radius);
	resultExact[5] = resultExact[1];

	for(size_t i = 0; i < resultExact.size(); ++i)
		ASSERT_NEAR(result[i], resultExact[i],
			std::numeric_limits<scalar_t>::epsilon());
}

// Sphere intersects one vertex (and thus 3 edges and 3 faces).
TEST(SphereTetAreaTest, Vertex) {
	vector_t v0{-std::sqrt(3) / 6.0, -1.0 / 2.0, 0};
	vector_t v1{std::sqrt(3) / 3.0, 0, 0};
	vector_t v2{-std::sqrt(3) / 6.0, +1.0 / 2.0, 0};
	vector_t v3{0, 0, std::sqrt(6) / 3.0};

	Tetrahedron tet{v0, v1, v2, v3};
	Sphere s({0, 0, 1.5}, 1.25);

	auto result = overlapArea(s, tet);

	// Compare with approximate result obtained via an Monte Carlo approach.
	std::array<scalar_t, 6> resultApprox;

	resultApprox[0] = scalar_t(0.19005658402860406);
	resultApprox[1] = scalar_t(0);
	resultApprox[2] = scalar_t(0.1879051823986737);
	resultApprox[3] = resultApprox[2];
	resultApprox[4] = resultApprox[2];
	resultApprox.back() = std::accumulate(resultApprox.begin() + 1,
		resultApprox.end() - 1, scalar_t(0));

	// Should be 1 / sqrt(N_{samples}), but does not quite work...
	constexpr scalar_t epsilon = 3.5e-06;
	for(size_t i = 0; i < resultApprox.size(); ++i)
		ASSERT_NEAR(result[i], resultApprox[i], epsilon);
}
