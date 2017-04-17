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

#define BOOST_TEST_MODULE sphere_element_overlap
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "overlap.hpp"

#include "common.hpp"

// Sphere intersects one edge (and thus 2 faces).
BOOST_AUTO_TEST_CASE(sphere_element_overlap_edge) {
	Sphere s(vector_t{0, -1, 1}, 1);

	scalar_t delta(std::sqrt(std::numeric_limits<scalar_t>::epsilon()));

	overlapWorker(s, unitHexahedron(), delta, 0.25 * s.volume);
}

// Sphere intersects one vertex (and thus 3 edges and 3 faces)
BOOST_AUTO_TEST_CASE(sphere_element_overlap_vertex) {
	Sphere s(vector_t{1, -1, 1}, 1);

	scalar_t delta(std::sqrt(std::numeric_limits<scalar_t>::epsilon()));

	overlapWorker(s, unitHexahedron(), delta, 0.125 * s.volume);
}

// Sphere outside of hexahedron, touching one face.
BOOST_AUTO_TEST_CASE(sphere_element_overlap_face_touching) {
	Sphere s(vector_t{0, 2, 0}, 1);

	scalar_t delta(std::sqrt(std::numeric_limits<scalar_t>::epsilon()));

	overlapWorker(s, unitHexahedron(), delta, 0);
}

// Sphere outside of hexahedron, touching one vertex.
BOOST_AUTO_TEST_CASE(sphere_element_overlap_vertex_touching) {
	Sphere s(vector_t{2, -1, 1}, 1);

	scalar_t delta(std::sqrt(std::numeric_limits<scalar_t>::epsilon()));

	overlapWorker(s, unitHexahedron(), delta, 0);
}

// Sphere outside of hexahedron, slightly overlapping one vertex.
BOOST_AUTO_TEST_CASE(sphere_element_overlap_vertex_edgecase) {
	Sphere s(vector_t{2 - 10 * detail::tinyEpsilon, -1, 1}, 1);

	scalar_t delta(std::sqrt(std::numeric_limits<scalar_t>::epsilon()));

	overlapWorker(s, unitHexahedron(), delta, 0);
}
