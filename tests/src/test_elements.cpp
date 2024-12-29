// Copyright (C) 2021-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

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
