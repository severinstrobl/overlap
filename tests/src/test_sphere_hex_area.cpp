// Copyright (C) 2015-2025 Severin Strobl <git@severin-strobl.de>
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
    const auto sphere = Sphere{{0, 0, 1}, 0.75};

    const auto result = overlap_area(sphere, hex);

    auto result_exact = std::array<Scalar, 8>{};
    result_exact[0] = Scalar{0.5} * sphere.surface_area();
    result_exact[6] = sphere.disk_area(sphere.radius);
    result_exact[7] = result_exact[6];

    for (auto i = 0u; i < result_exact.size(); ++i) {
      CHECK(result[i] == Approx(result_exact[i])
                             .epsilon(std::numeric_limits<Scalar>::epsilon()));
    }
  }

  // Sphere intersects one edge (and thus 1 edge and 2 faces).
  TEST_CASE("Edge") {
    const auto hex = unit_hexahedron();
    const auto sphere = Sphere{{1, 1, 0}, 0.75};

    const auto result = overlap_area(sphere, hex);

    auto result_exact = std::array<Scalar, 8>{};
    result_exact[0] = Scalar{0.25} * sphere.surface_area();
    result_exact[3] = 0.5 * sphere.disk_area(sphere.radius);
    result_exact[4] = result_exact[3];
    result_exact[7] = 2 * result_exact[3];

#if defined(__aarch64__)
    // A larger epsilon is required here as we seem to hit an edge case where
    // the calculation of the intersection points between the sphere and the
    // edge suffers from numerical inaccuracies when using GCC or Clang on
    // AArch64. For Clang, using either ffp-contract=fast or ffp-contract=on the
    // expected result is obtained.
    const auto epsilon = std::sqrt(std::numeric_limits<Scalar>::epsilon() *
                                   sphere.surface_area());
#else
    const auto epsilon = std::numeric_limits<Scalar>::epsilon();
#endif

    for (auto i = 0u; i < result_exact.size(); ++i) {
      INFO("result: ", result[i], ", expected: ", result_exact[i],
           ", delta:", result[i] - result_exact[i]);
      CHECK(result[i] == Approx(result_exact[i]).epsilon(epsilon));
    }
  }

  // Sphere intersects one edge (and thus 1 edge and 2 faces).
  TEST_CASE("EdgeOffCenter") {
    const auto hex = unit_hexahedron();
    const auto sphere = Sphere{{1.25, 0, 1}, 0.75};

    validate_overlap_area(
        sphere, hex,
        std::numeric_limits<Scalar>::epsilon() * sphere.surface_area());
  }

  // Sphere intersects one vertex (and thus 3 edge and 3 faces).
  TEST_CASE("Vertex") {
    const auto hex = unit_hexahedron();
    const auto sphere = Sphere{{1, 1, 1}, 0.75};

    const auto result = overlap_area(sphere, hex);

    auto result_exact = std::array<Scalar, 8>{};
    result_exact[0] = Scalar{0.125} * sphere.surface_area();
    result_exact[3] = 0.25 * sphere.disk_area(sphere.radius);
    result_exact[4] = result_exact[3];
    result_exact[6] = result_exact[3];
    result_exact[7] = 3 * result_exact[3];

    constexpr auto epsilon =
        Scalar{1e3} * std::numeric_limits<Scalar>::epsilon();

    for (auto i = 0u; i < result_exact.size(); ++i) {
      CHECK(result[i] == Approx(result_exact[i]).epsilon(epsilon));
    }
  }
}
