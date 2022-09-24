/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2021-2022 Severin Strobl <git@severin-strobl.de>
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

#include "common.hpp"

TEST_SUITE("contains") {
  using namespace overlap;

  TEST_CASE("SpherePoint") {
    CHECK(detail::contains(Sphere{}, Vector::Zero()));
    CHECK(detail::contains(Sphere{}, Vector::Constant(0.25)));
    CHECK(!detail::contains(Sphere{}, Vector::Constant(2.0)));
  }

  TEST_CASE("PolygonPoint") {
    const auto tri = detail::Triangle{{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}}};

    CHECK(detail::contains(tri, Vector{0.25, 0.25, 0.0}));
    CHECK(!detail::contains(tri, Vector{1.0, 1.0, 0.0}));

    const auto quad =
        detail::Quadrilateral{{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}}};

    CHECK(detail::contains(quad, Vector{0.5, 0.5, 0.0}));
    CHECK(!detail::contains(quad, Vector{-1.0, 0.0, 0.0}));
  }

  TEST_CASE("ElementPoint") {
    CHECK(detail::contains(unit_hexahedron(), Vector::Zero()));
    CHECK(detail::contains(unit_hexahedron(), Vector::Constant(0.5)));
    CHECK(!detail::contains(unit_hexahedron(), Vector::Constant(2.0)));
  }

  TEST_CASE("SphereElement") {
    CHECK(detail::contains(Sphere{Vector::Zero(), 3.0}, unit_hexahedron()));
    CHECK(!detail::contains(Sphere{Vector::Zero(), 0.5}, unit_hexahedron()));
  }
}
