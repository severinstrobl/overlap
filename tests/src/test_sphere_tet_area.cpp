// Copyright (C) 2015-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("SphereTetAreaTest") {
  using namespace overlap;

  // Sphere inside of element.
  TEST_CASE("SphereInTet") {
    const auto v0 = Vector{-std::sqrt(3) / 6.0, -1.0 / 2.0, 0};
    const auto v1 = Vector{std::sqrt(3) / 3.0, 0, 0};
    const auto v2 = Vector{-std::sqrt(3) / 6.0, +1.0 / 2.0, 0};
    const auto v3 = Vector{0, 0, std::sqrt(6) / 3.0};

    const auto tet = Tetrahedron{v0, v1, v2, v3};
    const auto s = Sphere({0, 0, 0.25}, 0.125);

    const auto result = overlap_area(s, tet);

    std::array<Scalar, 6> resultExact{};
    resultExact[0] = s.surface_area();

    for (std::size_t i = 0; i < resultExact.size(); ++i) {
      CHECK(result[i] == resultExact[i]);
    }
  }

  // Element contained in sphere.
  TEST_CASE("TetInSphere") {
    const auto v0 = Vector{-std::sqrt(3) / 6.0, -1.0 / 2.0, 0};
    const auto v1 = Vector{std::sqrt(3) / 3.0, 0, 0};
    const auto v2 = Vector{-std::sqrt(3) / 6.0, +1.0 / 2.0, 0};
    const auto v3 = Vector{0, 0, std::sqrt(6) / 3.0};

    const auto tet = Tetrahedron{v0, v1, v2, v3};
    const auto s = Sphere({0, 0, std::sqrt(6) / 6.0}, 2);

    const auto result = overlap_area(s, tet);

    std::array<Scalar, 6> resultExact{};
    resultExact[0] = Scalar{0};
    resultExact[1] = tet.faces[0].area;
    resultExact[2] = tet.faces[1].area;
    resultExact[3] = tet.faces[2].area;
    resultExact[4] = tet.faces[3].area;
    resultExact.back() = std::accumulate(resultExact.begin() + 1,
                                         resultExact.end() - 1, Scalar{0});

    for (std::size_t i = 0; i < resultExact.size(); ++i) {
      CHECK(result[i] == Approx(resultExact[i])
                             .epsilon(std::numeric_limits<Scalar>::epsilon()));
    }
  }

  // Sphere intersects one face.
  TEST_CASE("Face") {
    const auto v0 = Vector{-std::sqrt(3) / 6.0, -1.0 / 2.0, 0};
    const auto v1 = Vector{std::sqrt(3) / 3.0, 0, 0};
    const auto v2 = Vector{-std::sqrt(3) / 6.0, +1.0 / 2.0, 0};
    const auto v3 = Vector{0, 0, std::sqrt(6) / 3.0};

    const auto tet = Tetrahedron{v0, v1, v2, v3};
    const auto s = Sphere{{0, 0, 0}, 0.25};

    const auto result = overlap_area(s, tet);

    std::array<Scalar, 6> resultExact{};
    resultExact.fill(Scalar{0});
    resultExact[0] = Scalar{0.5} * s.surface_area();
    resultExact[1] = s.disk_area(s.radius);
    resultExact[5] = resultExact[1];

    for (std::size_t i = 0; i < resultExact.size(); ++i) {
      CHECK(result[i] == Approx(resultExact[i])
                             .epsilon(std::numeric_limits<Scalar>::epsilon()));
    }
  }

  // Sphere intersects one vertex (and thus 3 edges and 3 faces).
  TEST_CASE("Vertex") {
    const auto v0 = Vector{-std::sqrt(3) / 6.0, -1.0 / 2.0, 0};
    const auto v1 = Vector{std::sqrt(3) / 3.0, 0, 0};
    const auto v2 = Vector{-std::sqrt(3) / 6.0, +1.0 / 2.0, 0};
    const auto v3 = Vector{0, 0, std::sqrt(6) / 3.0};

    const auto tet = Tetrahedron{v0, v1, v2, v3};
    const auto s = Sphere{{0, 0, 1.5}, 1.25};

    const auto result = overlap_area(s, tet);

    // Compare with approximate result obtained via an Monte Carlo approach.
    std::array<Scalar, 6> resultApprox{};
    resultApprox[0] = Scalar{0.19005658402860406};
    resultApprox[1] = Scalar{0};
    resultApprox[2] = Scalar{0.1879051823986737};
    resultApprox[3] = resultApprox[2];
    resultApprox[4] = resultApprox[2];
    resultApprox.back() = std::accumulate(resultApprox.begin() + 1,
                                          resultApprox.end() - 1, Scalar{0});

    // Should be 1 / sqrt(N_{samples}), but does not quite work...
    constexpr auto epsilon = Scalar{3.5e-06};
    for (std::size_t i = 0; i < resultApprox.size(); ++i) {
      CHECK(result[i] == Approx(resultApprox[i]).epsilon(epsilon));
    }
  }
}
