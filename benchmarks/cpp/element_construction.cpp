// Copyright (C) 2026 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include <array>
#include <cmath>

#include <nanobench.h>

#include <doctest/doctest.h>

#include "overlap/overlap.hpp"

#include "common.hpp"

TEST_SUITE("ElementConstruction") {
  using overlap::Hexahedron;
  using overlap::Scalar;
  using overlap::Tetrahedron;
  using overlap::Vector;
  using overlap::Wedge;

  TEST_CASE("Tetrahedron") {
    const auto sqrt3 = std::sqrt(Scalar{3});
    const auto sqrt6 = std::sqrt(Scalar{6});

    // clang-format off
    const auto tet_vertices = std::array<Vector, 4>{{
      {-sqrt3 / 6.0, -0.5, 0}, {sqrt3 / 3.0, 1.0, 0},
      {-sqrt3 / 6.0,  0.5, 0}, {0, 0, sqrt6 / 3.0}}};
    // clang-format on

    create_benchmark("construction[tet]", [&tet_vertices]() {
      const auto result = Tetrahedron{tet_vertices};
      ankerl::nanobench::doNotOptimizeAway(result);
    });
  }

  TEST_CASE("Wedge") {
    // clang-format off
    const auto wedge_vertices = std::array<Vector, 6>{{
      {0, 0, -1}, {1, 0, -1}, {0, 1, -1},
      {0, 0,  1}, {1, 0,  1}, {0, 1,  1}}};
    // clang-format on

    create_benchmark("construction[wedge]", [&wedge_vertices]() {
      const auto result = Wedge{wedge_vertices};
      ankerl::nanobench::doNotOptimizeAway(result);
    });
  }

  TEST_CASE("Hexahedron") {
    // clang-format off
    const auto hex_vertices = std::array<Vector, 8>{{
      {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
      {-1, -1,  1}, {1, -1,  1}, {1, 1,  1}, {-1, 1,  1}}};
    // clang-format on

    create_benchmark("construction[hex]", [&hex_vertices]() {
      const auto result = Hexahedron{hex_vertices};
      ankerl::nanobench::doNotOptimizeAway(result);
    });
  }
}
