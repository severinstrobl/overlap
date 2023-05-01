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

TEST_SUITE("Elements") {
  using namespace overlap;
  constexpr auto epsilon = std::numeric_limits<Scalar>::epsilon();

  TEST_CASE("Tetrahedron") {
    SUBCASE("DefaultConstructor") {
      const auto tet0 = Tetrahedron{};
      REQUIRE(tet0.volume == Scalar{0});
    }

    const auto sqrt3 = std::sqrt(Scalar{3});
    const auto sqrt6 = std::sqrt(Scalar{6});

    SUBCASE("ConstructFromInitializerList") {
      // clang-format off
      const auto tet = Tetrahedron{{{
        {-sqrt3 / 6.0, -1.0 / 2.0, 0}, {sqrt3 / 3.0, 0, 0},
        {-sqrt3 / 6.0, +1.0 / 2.0, 0}, {0, 0, sqrt6 / 3.0}}}};
      // clang-format on

      CHECK(tet.volume ==
            Approx((1.0 / 12.0) * std::sqrt(2.0)).epsilon(epsilon));
    }

    // clang-format off
    const auto vertices = std::array<Vector, 4>{{
      {-sqrt3 / 6.0, -1.0 / 2.0, 0}, {sqrt3 / 3.0, 0, 0},
      {-sqrt3 / 6.0, +1.0 / 2.0, 0}, {0, 0, sqrt6 / 3.0}}};
    // clang-format on

    SUBCASE("ConstructFromArray") {
      CHECK(Tetrahedron{vertices}.volume ==
            Approx((1.0 / 12.0) * std::sqrt(2.0)).epsilon(epsilon));
    }

    SUBCASE("InvalidNodeOrder") {
      REQUIRE_THROWS_AS(
          Tetrahedron(vertices[0], vertices[1], vertices[3], vertices[2]),
          AssertionError);
    }
  }

  TEST_CASE("Wedge") {
    SUBCASE("DefaultConstructor") {
      const auto wedge = Wedge{};
      REQUIRE(wedge.volume == Scalar{0});
    }

    SUBCASE("ConstructFromInitializerList") {
      // clang-format off
      const auto wedge = Wedge{{{
        {-1, -1, -1}, {1, -1, -1}, {1, 1, -1},
        {-1, -1, 1},  {1, -1, 1},  {1, 1, 1}}}};
      // clang-format on

      CHECK(wedge.volume == Approx(Scalar{4}).epsilon(epsilon));
    }

    SUBCASE("ConstructFromArray") {
      // clang-format off
      const auto wedge = Wedge{std::array<Vector, 6>{{
        {-1, -1, -1}, {1, -1, -1}, {1, 1, -1},
        {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}}}};
      // clang-format on

      CHECK(wedge.volume == Approx(Scalar{4}).epsilon(epsilon));
    }
  }

  TEST_CASE("Hexahedron") {
    SUBCASE("DefaultConstructor") {
      const auto hex = Hexahedron{};
      REQUIRE(hex.volume == Scalar{0});
    }

    SUBCASE("ConstructFromInitializerList") {
      // clang-format off
      const auto hex = Hexahedron{{{
        {-1, -1, -1}, {1, -1, -1},
        { 1,  1, -1}, {-1,  1, -1},

        {-1, -1,  1}, {1, -1,  1},
        { 1,  1,  1}, {-1,  1,  1}}}};
      // clang-format on

      CHECK(hex.volume == Approx(Scalar{8}).epsilon(epsilon));
    }

    SUBCASE("ConstructFromArray") {
      // clang-format off

      // clang-format off
      const auto hex = Hexahedron{std::array<Vector, 8>{{
        {-1, -1, -1}, {1, -1, -1}, {1,  1, -1}, {-1,  1, -1},
        {-1, -1,  1}, {1, -1,  1}, {1,  1,  1}, {-1,  1,  1}}}};
      // clang-format on

      CHECK(hex.volume == Approx(Scalar{8}).epsilon(epsilon));
    }
  }
}