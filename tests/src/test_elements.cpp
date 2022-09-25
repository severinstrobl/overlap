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
  TEST_CASE("Tetrahedron") {
    using namespace overlap;

    const auto tet0 = Tetrahedron{};
    REQUIRE(tet0.volume == Scalar{0});

    const auto sqrt3 = std::sqrt(Scalar{3});
    const auto sqrt6 = std::sqrt(Scalar{6});

    // clang-format off
    const auto tet1 = Tetrahedron{{{
      {-sqrt3 / 6.0, -1.0 / 2.0, 0}, {sqrt3 / 3.0, 0, 0},
      {-sqrt3 / 6.0, +1.0 / 2.0, 0}, {0, 0, sqrt6 / 3.0}}}};

    const auto tet2 = Tetrahedron{std::array<Vector, 4>{{
      {-sqrt3 / 6.0, -1.0 / 2.0, 0}, {sqrt3 / 3.0, 0, 0},
      {-sqrt3 / 6.0, +1.0 / 2.0, 0}, {0, 0, sqrt6 / 3.0}}}};
    // clang-format on

    constexpr auto epsilon = std::numeric_limits<Scalar>::epsilon();

    CHECK(tet1.volume ==
          Approx((1.0 / 12.0) * std::sqrt(2.0)).epsilon(epsilon));
    CHECK(tet1.volume == Approx(tet2.volume).epsilon(epsilon));
  }

  TEST_CASE("Wedge") {
    using namespace overlap;

    const auto wedge0 = Wedge{};
    REQUIRE(wedge0.volume == Scalar{0});

    // clang-format off
    const auto wedge1 = Wedge{{{
      {-1, -1, -1}, {1, -1, -1}, {1, 1, -1},
      {-1, -1, 1},  {1, -1, 1},  {1, 1, 1}}}};

    const auto wedge2 = Wedge{std::array<Vector, 6>{{
      {-1, -1, -1}, {1, -1, -1}, {1, 1, -1},
      {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}}}};
    // clang-format on

    constexpr auto epsilon = std::numeric_limits<Scalar>::epsilon();
    CHECK(wedge1.volume == Approx(Scalar{4}).epsilon(epsilon));
    CHECK(wedge1.volume == Approx(wedge2.volume).epsilon(epsilon));
  }

  TEST_CASE("Hexahedron") {
    using namespace overlap;

    const auto hex0 = Hexahedron{};
    REQUIRE(hex0.volume == Scalar{0});

    // clang-format off
    const auto hex1 = Hexahedron{{{
      {-1, -1, -1}, {1, -1, -1},
      { 1,  1, -1}, {-1,  1, -1},

      {-1, -1,  1}, {1, -1,  1},
      { 1,  1,  1}, {-1,  1,  1}}}};

    const auto hex2 = Hexahedron{std::array<Vector, 8>{{
      {-1, -1, -1}, {1, -1, -1}, {1,  1, -1}, {-1,  1, -1},
      {-1, -1,  1}, {1, -1,  1}, {1,  1,  1}, {-1,  1,  1}}}};
    // clang-format on

    constexpr auto epsilon = std::numeric_limits<Scalar>::epsilon();
    CHECK(hex1.volume == Approx(Scalar{8}).epsilon(epsilon));
    CHECK(hex1.volume == Approx(hex2.volume).epsilon(epsilon));
  }
}