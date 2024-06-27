// Copyright (C) 2015-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("DecomposeElements") {
  TEST_CASE("Hexahedron") {
    using namespace overlap;

    const auto hex = overlap::unit_hexahedron();
    std::array<Tetrahedron, 5> tets5;
    std::array<Tetrahedron, 6> tets6;
    std::array<Wedge, 2> wedges;

    decompose(hex, tets5);
    decompose(hex, tets6);
    decompose(hex, wedges);

    auto tets5Volume = Scalar{0};
    for (const auto& tet : tets5) {
      tets5Volume += tet.volume;
    }

    auto tets6Volume = Scalar{0};
    auto tets24Volume = Scalar{0};
    for (const auto& tet : tets6) {
      std::array<Tetrahedron, 4> subTets;
      decompose(tet, subTets);
      tets6Volume += tet.volume;

      for (const auto& subTet : subTets) {
        tets24Volume += subTet.volume;
      }
    }

    constexpr auto epsilon =
        Scalar{5e2} * std::numeric_limits<Scalar>::epsilon();

    REQUIRE(hex.volume == Approx(tets5Volume).epsilon(epsilon));
    REQUIRE(hex.volume == Approx(tets6Volume).epsilon(epsilon));
    REQUIRE(hex.volume == Approx(tets24Volume).epsilon(epsilon));
    REQUIRE(hex.volume ==
            Approx(wedges[0].volume + wedges[1].volume).epsilon(epsilon));
  }
}
