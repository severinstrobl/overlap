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

static const scalar_t epsilon =
	std::sqrt(std::numeric_limits<scalar_t>::epsilon());

// Sphere outside of hexahedron, touching one face.
TEST(SphereElementOverlap, Face) {
	Sphere s(vector_t{0, 2, 0}, 1);

	overlap(s, unitHexahedron(), epsilon, scalar_t{0});
}

// Sphere intersects one edge (and thus 2 faces).
TEST(SphereElementOverlap, Edge) {
	Sphere s(vector_t{0, -1, 1}, 1);

	overlap(s, unitHexahedron(), epsilon, 0.25 * s.volume);
}

// Sphere intersects one vertex (and thus 3 edges and 3 faces)
TEST(SphereElementOverlap, Vertex) {
	Sphere s(vector_t{1, -1, 1}, 1);

	overlap(s, unitHexahedron(), epsilon, 0.125 * s.volume);
}

// Sphere outside of hexahedron, touching one vertex.
TEST(SphereElementOverlap, VertexTouching) {
	Sphere s(vector_t{2, -1, 1}, 1);

	overlap(s, unitHexahedron(), epsilon, scalar_t{0});
}

// Sphere outside of hexahedron, slightly overlapping one vertex.
TEST(SphereElementOverlap, Vertex2) {
	Sphere s(vector_t{2 - 10 * detail::tinyEpsilon, -1, 1}, 1);

	overlap(s, unitHexahedron(), epsilon, scalar_t{0});
}
