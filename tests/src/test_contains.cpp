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

#include "gtest/gtest.h"

#include "overlap/overlap.hpp"

#include "common.hpp"

TEST(Contains, SpherePoint) {
  using namespace overlap;

  ASSERT_TRUE(detail::contains(Sphere{}, Vector::Zero()));
  ASSERT_TRUE(detail::contains(Sphere{}, Vector::Constant(0.25)));
  ASSERT_FALSE(detail::contains(Sphere{}, Vector::Constant(2.0)));
}

TEST(Contains, PolygonPoint) {
  using namespace overlap;

  const auto tri = detail::Triangle{{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}}};

  ASSERT_TRUE(detail::contains(tri, Vector{0.25, 0.25, 0.0}));
  ASSERT_FALSE(detail::contains(tri, Vector{1.0, 1.0, 0.0}));

  const auto quad =
      detail::Quadrilateral{{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}}};

  ASSERT_TRUE(detail::contains(quad, Vector{0.5, 0.5, 0.0}));
  ASSERT_FALSE(detail::contains(quad, Vector{-1.0, 0.0, 0.0}));
}

TEST(Contains, ElementPoint) {
  using namespace overlap;

  ASSERT_TRUE(detail::contains(unit_hexahedron(), Vector::Zero()));
  ASSERT_TRUE(detail::contains(unit_hexahedron(), Vector::Constant(0.5)));
  ASSERT_FALSE(detail::contains(unit_hexahedron(), Vector::Constant(2.0)));
}

TEST(Contains, SphereElement) {
  using namespace overlap;

  ASSERT_TRUE(detail::contains(Sphere{Vector::Zero(), 3.0}, unit_hexahedron()));
  ASSERT_FALSE(detail::contains(Sphere{Vector::Zero(), 0.5}, unit_hexahedron()));
}