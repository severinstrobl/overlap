/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2015-2022 Severin Strobl <severin.strobl@fau.de>
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

#include "overlap/overlap.hpp"

#include "common.hpp"

static const auto epsilon = std::sqrt(std::numeric_limits<Scalar>::epsilon());

// Sphere outside of hexahedron, touching one face.
TEST(SphereElementOverlap, Face) {
  const auto sphere = Sphere{{0, 2, 0}, 1};

  validate_overlap_volume(sphere, unit_hexahedron(), epsilon, Scalar{0});
}

// Sphere intersects one edge (and thus 2 faces).
TEST(SphereElementOverlap, Edge) {
  const auto sphere = Sphere{Vector{0, -1, 1}, 1};

  validate_overlap_volume(sphere, unit_hexahedron(), epsilon,
                          0.25 * sphere.volume);
}

// Sphere intersects one edge (and thus 2 faces), edge passing through center
// of sphere -> spherical wedge with angle pi/4.
TEST(SphereElementOverlap, Wedge) {
  const auto sphere = Sphere{Vector::Zero(), 1};

  // clang-format off
  const auto hex = Hexahedron{{{
    {0, 0, -1}, {2, 2, -1},
    {2, 4, -1}, {0, 4, -1},

    {0, 0, 1}, {2, 2, 1},
    {2, 4, 1}, {0, 4, 1}}}};
  // clang-format on

  validate_overlap_volume(
      sphere, hex, epsilon,
      2.0 / 3.0 * sphere.radius * sphere.radius * 0.25 * detail::pi);
}

// Sphere intersects one vertex (and thus 3 edges and 3 faces)
TEST(SphereElementOverlap, Vertex) {
  const auto sphere = Sphere{Vector{1, -1, 1}, 1};

  validate_overlap_volume(sphere, unit_hexahedron(), epsilon,
                          0.125 * sphere.volume);
}

// Sphere outside of hexahedron, touching one vertex.
TEST(SphereElementOverlap, VertexTouching) {
  const auto sphere = Sphere{Vector{2, -1, 1}, 1};

  validate_overlap_volume(sphere, unit_hexahedron(), epsilon, Scalar{0});
}

// Sphere outside of hexahedron, slightly overlapping one vertex.
TEST(SphereElementOverlap, VertexOverlap) {
  const auto sphere = Sphere{Vector{2 - 10 * detail::tinyEpsilon, -1, 1}, 1};

  validate_overlap_volume(sphere, unit_hexahedron(), epsilon, Scalar{0});
}
