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

#include <doctest/doctest.h>

#include <array>
#include <numeric>

#include "common.hpp"
#include "overlap/overlap.hpp"

namespace overlap {

TEST_SUITE("Polyhedron") {
  using namespace detail;

  const auto sqrt3 = std::sqrt(Scalar{3});
  const auto sqrt6 = std::sqrt(Scalar{6});

  TEST_CASE("ConstructorInitializerList") {
    const auto poly = Tetrahedron{{{{-sqrt3 / 6.0, -1.0 / 2.0, 0},
                                    {sqrt3 / 3.0, 0, 0},
                                    {-sqrt3 / 6.0, 1.0 / 2.0, 0},
                                    {0, 0, sqrt6 / 3.0}}}};

    CHECK_EQ(poly.vertices(),
             std::array{Vector{-sqrt3 / 6, -0.5, 0}, Vector{sqrt3 / 3, 0, 0},
                        Vector{-sqrt3 / 6, 0.5, 0}, Vector{0, 0, sqrt6 / 3}});
  }

  TEST_CASE("ConstructorInitializerListMixedTypes") {
    const auto poly = Polygon<3>{{{{0.0, 0, 0}, {1, 0, 0}, {1.0, 1, 0.0}}}};

    CHECK(poly.vertices() ==
          std::array{Vector{0, 0, 0}, Vector{1, 0, 0}, Vector{1, 1, 0}});
  }

  TEST_CASE("ConstructorArrayOfVectors") {
    const auto poly = Polygon<3>{
        std::array{Vector{0, 0, 0}, Vector{1, 0, 0}, Vector{1, 1, 0}}};

    CHECK(poly.vertices() ==
          std::array{Vector{0, 0, 0}, Vector{1, 0, 0}, Vector{1, 1, 0}});
  }

  constexpr auto epsilon = std::numeric_limits<Scalar>::epsilon();

  TEST_CASE("Tetrahedron") {
    const auto vertices = std::array{
        Vector{-sqrt3 / 6.0, -1.0 / 2.0, 0}, Vector{sqrt3 / 3.0, 0, 0},
        Vector{-sqrt3 / 6.0, +1.0 / 2.0, 0}, Vector{0, 0, sqrt6 / 3.0}};

    SUBCASE("InvalidVertexOrdering") {
      CHECK_THROWS_WITH_AS(
          Tetrahedron(vertices[0], vertices[1], vertices[2], -vertices[3]),
          doctest::Contains("invalid vertex order detected"), AssertionError);
    }

    SUBCASE("Center") {
      const auto tet = Tetrahedron{vertices};
      const auto center =
          (1.0 / 4.0 *
           std::accumulate(std::begin(vertices), std::end(vertices),
                           Vector::Zero().eval()))
              .eval();

      CHECK_EQ(tet.center(), center);
    }

    SUBCASE("Volume") {
      const auto tet = Tetrahedron{vertices};

      CHECK_EQ(tet.volume(),
               Approx((1.0 / 12.0) * std::sqrt(2.0)).epsilon(epsilon));
    }

    SUBCASE("Faces") {
      const auto tet = Tetrahedron{vertices};

      // face #0
      CHECK_EQ(tet.face(0).vertex(0), vertices[0]);
      CHECK_EQ(tet.face(0).vertex(1), vertices[2]);
      CHECK_EQ(tet.face(0).vertex(2), vertices[1]);

      CHECK_EQ(tet.face(0).normal(), -Vector::UnitZ());

      // face #1
      CHECK_EQ(tet.face(1).vertex(0), vertices[0]);
      CHECK_EQ(tet.face(1).vertex(1), vertices[1]);
      CHECK_EQ(tet.face(1).vertex(2), vertices[3]);

      // face #2
      CHECK_EQ(tet.face(2).vertex(0), vertices[1]);
      CHECK_EQ(tet.face(2).vertex(1), vertices[2]);
      CHECK_EQ(tet.face(2).vertex(2), vertices[3]);

      // face #3
      CHECK_EQ(tet.face(3).vertex(0), vertices[2]);
      CHECK_EQ(tet.face(3).vertex(1), vertices[0]);
      CHECK_EQ(tet.face(3).vertex(2), vertices[3]);
    }

    SUBCASE("Decompose") {
      const auto tet = Tetrahedron{vertices};
      auto sub_tets = std::array<Tetrahedron, 4>{};
      decompose(tet, sub_tets);

      CHECK_EQ(Approx(tet.volume()),
               std::accumulate(std::begin(sub_tets), std::end(sub_tets),
                               Scalar{0}, [](const auto sum, const auto& tet) {
                                 return sum + tet.volume();
                               }));
    }
  }

  TEST_CASE("Wedge") {
    const auto vertices =
        std::array{Vector{-1, -1, -1}, Vector{1, -1, -1}, Vector{1, 1, -1},
                   Vector{-1, -1, 1},  Vector{1, -1, 1},  Vector{1, 1, 1}};

    SUBCASE("Center") {
      const auto wedge = Wedge{vertices};
      const auto center =
          (1.0 / 6.0 *
           std::accumulate(std::begin(vertices), std::end(vertices),
                           Vector::Zero().eval()))
              .eval();

      CHECK_EQ(wedge.center(), center);
    }

    SUBCASE("Volume") {
      const auto wedge = Wedge{vertices};

      CHECK_EQ(wedge.volume(), Approx(4.0).epsilon(epsilon));
    }

    SUBCASE("Faces") {
      const auto wedge = Wedge{vertices};

      // face #0
      CHECK_EQ(wedge.face(0).vertex(0), vertices[0]);
      CHECK_EQ(wedge.face(0).vertex(1), vertices[1]);
      CHECK_EQ(wedge.face(0).vertex(2), vertices[4]);
      CHECK_EQ(wedge.face(0).vertex(3), vertices[3]);

      CHECK_EQ(wedge.face(0).normal(), -Vector::UnitY());

      // face #1
      CHECK_EQ(wedge.face(1).vertex(0), vertices[1]);
      CHECK_EQ(wedge.face(1).vertex(1), vertices[2]);
      CHECK_EQ(wedge.face(1).vertex(2), vertices[5]);
      CHECK_EQ(wedge.face(1).vertex(3), vertices[4]);

      CHECK_EQ(wedge.face(1).normal(), Vector::UnitX());

      // face #2
      CHECK_EQ(wedge.face(2).vertex(0), vertices[2]);
      CHECK_EQ(wedge.face(2).vertex(1), vertices[0]);
      CHECK_EQ(wedge.face(2).vertex(2), vertices[3]);
      CHECK_EQ(wedge.face(2).vertex(3), vertices[5]);

      CHECK_EQ(wedge.face(2).normal(),
               (-Vector::UnitX() + Vector::UnitY()).normalized());

      // face #3
      CHECK_EQ(wedge.face(3).vertex(0), vertices[0]);
      CHECK_EQ(wedge.face(3).vertex(1), vertices[2]);
      CHECK_EQ(wedge.face(3).vertex(2), vertices[1]);

      CHECK_EQ(wedge.face(3).normal(), -Vector::UnitZ());

      // face #4
      CHECK_EQ(wedge.face(4).vertex(0), vertices[3]);
      CHECK_EQ(wedge.face(4).vertex(1), vertices[4]);
      CHECK_EQ(wedge.face(4).vertex(2), vertices[5]);

      CHECK_EQ(wedge.face(4).normal(), Vector::UnitZ());
    }
  }

  TEST_CASE("Hexahedron") {
    const auto vertices =
        std::array{Vector{-1, -1, -1}, Vector{1, -1, -1}, Vector{1, 1, -1},
                   Vector{-1, 1, -1},  Vector{-1, -1, 1}, Vector{1, -1, 1},
                   Vector{1, 1, 1},    Vector{-1, 1, 1}};

    SUBCASE("Center") {
      const auto hex = Hexahedron{vertices};
      const auto center =
          (1.0 / 6.0 *
           std::accumulate(std::begin(vertices), std::end(vertices),
                           Vector::Zero().eval()))
              .eval();

      CHECK_EQ(hex.center(), center);
    }

    SUBCASE("Volume") {
      const auto hex = Hexahedron{vertices};

      CHECK_EQ(hex.volume(), Approx(8.0).epsilon(epsilon));
    }

    SUBCASE("Faces") {
      const auto hex = Hexahedron{vertices};

      // face #0
      CHECK_EQ(hex.face(0).vertex(0), vertices[0]);
      CHECK_EQ(hex.face(0).vertex(1), vertices[3]);
      CHECK_EQ(hex.face(0).vertex(2), vertices[2]);
      CHECK_EQ(hex.face(0).vertex(3), vertices[1]);

      CHECK_EQ(hex.face(0).normal(), -Vector::UnitZ());

      // face #1
      CHECK_EQ(hex.face(1).vertex(0), vertices[0]);
      CHECK_EQ(hex.face(1).vertex(1), vertices[1]);
      CHECK_EQ(hex.face(1).vertex(2), vertices[5]);
      CHECK_EQ(hex.face(1).vertex(3), vertices[4]);

      CHECK_EQ(hex.face(1).normal(), -Vector::UnitY());

      // face #2
      CHECK_EQ(hex.face(2).vertex(0), vertices[1]);
      CHECK_EQ(hex.face(2).vertex(1), vertices[2]);
      CHECK_EQ(hex.face(2).vertex(2), vertices[6]);
      CHECK_EQ(hex.face(2).vertex(3), vertices[5]);

      CHECK_EQ(hex.face(2).normal(), Vector::UnitX());

      // face #3
      CHECK_EQ(hex.face(3).vertex(0), vertices[2]);
      CHECK_EQ(hex.face(3).vertex(1), vertices[3]);
      CHECK_EQ(hex.face(3).vertex(2), vertices[7]);
      CHECK_EQ(hex.face(3).vertex(3), vertices[6]);

      CHECK_EQ(hex.face(3).normal(), Vector::UnitY());

      // face #4
      CHECK_EQ(hex.face(4).vertex(0), vertices[0]);
      CHECK_EQ(hex.face(4).vertex(1), vertices[4]);
      CHECK_EQ(hex.face(4).vertex(2), vertices[7]);
      CHECK_EQ(hex.face(4).vertex(3), vertices[3]);

      CHECK_EQ(hex.face(4).normal(), -Vector::UnitX());

      // face #5
      CHECK_EQ(hex.face(5).vertex(0), vertices[4]);
      CHECK_EQ(hex.face(5).vertex(1), vertices[5]);
      CHECK_EQ(hex.face(5).vertex(2), vertices[6]);
      CHECK_EQ(hex.face(5).vertex(3), vertices[7]);

      CHECK_EQ(hex.face(5).normal(), Vector::UnitZ());
    }

    SUBCASE("DecomposeIntoWedges") {
      const auto hex = Hexahedron{vertices};
      auto wedges = std::array<Wedge, 2>{};
      decompose(hex, wedges);

      CHECK_EQ(Approx(hex.volume()),
               std::accumulate(std::begin(wedges), std::end(wedges), Scalar{0},
                               [](const auto sum, const auto& wedge) {
                                 return sum + wedge.volume();
                               }));
    }

    SUBCASE("DecomposeInto5Tets") {
      const auto hex = Hexahedron{vertices};
      auto tets = std::array<Tetrahedron, 5>{};
      decompose(hex, tets);

      CHECK_EQ(Approx(hex.volume()),
               std::accumulate(std::begin(tets), std::end(tets), Scalar{0},
                               [](const auto sum, const auto& tet) {
                                 return sum + tet.volume();
                               }));
    }

    SUBCASE("DecomposeInto6Tets") {
      const auto hex = Hexahedron{vertices};
      auto tets = std::array<Tetrahedron, 6>{};
      decompose(hex, tets);

      CHECK_EQ(Approx(hex.volume()),
               std::accumulate(std::begin(tets), std::end(tets), Scalar{0},
                               [](const auto sum, const auto& tet) {
                                 return sum + tet.volume();
                               }));
    }
  }
}

}  // namespace overlap
