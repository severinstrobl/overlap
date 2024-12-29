// Copyright (C) 2015-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("SphereHexAreaTest") {
  using namespace overlap;

  // Sphere intersects one face.
  TEST_CASE("Face") {
    const auto hex = unit_hexahedron();
    const auto s = Sphere{{0, 0, 1}, 0.75};

    auto result = overlap_area(s, hex);

    std::array<Scalar, 8> resultExact{};
    resultExact[0] = Scalar{0.5} * s.surface_area();
    resultExact[6] = s.disk_area(s.radius);
    resultExact[7] = resultExact[6];

    for (std::size_t i = 0; i < resultExact.size(); ++i) {
      CHECK(result[i] == Approx(resultExact[i])
                             .epsilon(std::numeric_limits<Scalar>::epsilon()));
    }
  }

  // Sphere intersects one edge (and thus 1 edge and 2 faces).
  TEST_CASE("Edge") {
    const auto hex = unit_hexahedron();
    const auto sphere = Sphere{{1, 1, 0}, 0.75};

    auto result = overlap_area(sphere, hex);

    std::array<Scalar, 8> resultExact{};
    resultExact[0] = Scalar{0.25} * sphere.surface_area();
    resultExact[3] = 0.5 * sphere.disk_area(sphere.radius);
    resultExact[4] = resultExact[3];
    resultExact[7] = 2 * resultExact[3];

    for (std::size_t i = 0; i < resultExact.size(); ++i) {
      CHECK(result[i] == Approx(resultExact[i])
                             .epsilon(std::numeric_limits<Scalar>::epsilon()));
    }
  }

  // Sphere intersects one vertex (and thus 3 edge and 3 faces).
  TEST_CASE("Vertex") {
    const auto hex = unit_hexahedron();
    Sphere s({1, 1, 1}, 0.75);

    auto result = overlap_area(s, hex);

    std::array<Scalar, 8> resultExact{};
    resultExact[0] = Scalar{0.125} * s.surface_area();
    resultExact[3] = 0.25 * s.disk_area(s.radius);
    resultExact[4] = resultExact[3];
    resultExact[6] = resultExact[3];
    resultExact[7] = 3 * resultExact[3];

    constexpr auto epsilon =
        Scalar{1e3} * std::numeric_limits<Scalar>::epsilon();

    for (std::size_t i = 0; i < resultExact.size(); ++i) {
      CHECK(result[i] == Approx(resultExact[i]).epsilon(epsilon));
    }
  }
}
