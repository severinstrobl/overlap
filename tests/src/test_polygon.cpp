/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2021-2024 Severin Strobl <git@severin-strobl.de>
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

#include "overlap/overlap.hpp"

#include <array>
#include <iterator>
#include <numeric>

namespace overlap {

TEST_SUITE("Polygon") {
  using namespace detail;

  TEST_CASE("ConstructorInitializerList") {
    const auto poly =
        Polygon<3>{{{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}}};

    CHECK_EQ(poly.vertices(),
             std::array{Vector{0, 0, 0}, Vector{1, 0, 0}, Vector{1, 1, 0}});
  }

  TEST_CASE("ConstructorInitializerListMixedTypes") {
    const auto poly = Polygon<3>{{{{0.0, 0, 0}, {1, 0, 0}, {1.0, 1, 0.0}}}};

    CHECK_EQ(poly.vertices(),
             std::array{Vector{0, 0, 0}, Vector{1, 0, 0}, Vector{1, 1, 0}});
  }

  TEST_CASE("ConstructorArrayOfVectors") {
    const auto poly = Polygon<3>{
        std::array{Vector{0, 0, 0}, Vector{1, 0, 0}, Vector{1, 1, 0}}};

    CHECK_EQ(poly.vertices(),
             std::array{Vector{0, 0, 0}, Vector{1, 0, 0}, Vector{1, 1, 0}});
  }
}

TEST_SUITE("Triangle") {
  using namespace detail;

  const auto vertices =
      std::array{Vector{0, 0, 0}, Vector{1, 0, 0}, Vector{1, 1, 0}};

  TEST_CASE("Center") {
    const auto tri = Triangle{vertices};
    const auto center =
        (1.0 / 3.0 *
         std::accumulate(std::begin(vertices), std::end(vertices),
                         Vector::Zero().eval()))
            .eval();

    CHECK_EQ(tri.center(), center);
  }

  TEST_CASE("Normal") {
    const auto tri = Triangle{vertices};

    CHECK_EQ(tri.normal(), Vector::UnitZ());
  }

  TEST_CASE("Area") {
    const auto tri = Triangle{vertices};

    CHECK_EQ(tri.area(), 0.5);
  }

  TEST_CASE("IsPlanar") {
    const auto tri = Triangle{vertices};

    CHECK(tri.is_planar());
  }
}

TEST_SUITE("Quadrilateral") {
  using namespace detail;

  const auto vertices = std::array{Vector{0, 0, 0}, Vector{1, 0, 0},
                                   Vector{1, 1, 0}, Vector{0, 1, 0}};

  TEST_CASE("Center") {
    const auto quad = Quadrilateral{vertices};
    const auto center =
        (1.0 / 4.0 *
         std::accumulate(std::begin(vertices), std::end(vertices),
                         Vector::Zero().eval()))
            .eval();

    CHECK_EQ(quad.center(), center);
  }

  TEST_CASE("Normal") {
    const auto quad = Quadrilateral{vertices};

    CHECK_EQ(quad.normal(), Vector::UnitZ());
  }

  TEST_CASE("Area") {
    const auto quad = Quadrilateral{vertices};

    CHECK_EQ(quad.area(), 1.0);
  }

  TEST_CASE("IsPlanar") {
    const auto planar_quad = Quadrilateral{vertices};

    CHECK(planar_quad.is_planar());

    const auto non_planar_quad = Quadrilateral{
        {{vertices[0], vertices[1], vertices[2], {0.0, 1.0, 0.1}}}};

    CHECK_FALSE(non_planar_quad.is_planar());
  }
}

}  // namespace overlap
