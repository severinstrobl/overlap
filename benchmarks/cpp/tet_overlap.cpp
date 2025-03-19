// Copyright (C) 2023 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include <doctest/doctest.h>

#include "overlap/overlap.hpp"

#include "common.hpp"

TEST_SUITE("TetOverlap") {
  using namespace overlap;

  const auto sqrt3 = std::sqrt(Scalar{3});
  const auto sqrt6 = std::sqrt(Scalar{6});

  // clang-format off
  const auto tet = Tetrahedron{{{
    {-sqrt3 / 6.0, -1.0 / 2.0, 0}, {sqrt3 / 3.0, 0, 0},
    {-sqrt3 / 6.0, +1.0 / 2.0, 0}, {0, 0, sqrt6 / 3.0}}}};
  // clang-format on

  TEST_CASE("TetOverlapVolume") {
    const auto sphere = Sphere{Vector::Zero(), 0.5};

    create_benchmark("tet_overlap_volume[sphere-in-hex]", [&]() {
      const auto result = overlap_volume(sphere, tet);
      ankerl::nanobench::doNotOptimizeAway(result);
    });
  }

  TEST_CASE("TetOverlapVolume") {
    const auto sphere = Sphere{Vector::Zero(), 5.0};

    create_benchmark("tet_overlap_volume[hex-in-sphere]", [&]() {
      const auto result = overlap_volume(sphere, tet);
      ankerl::nanobench::doNotOptimizeAway(result);
    });
  }

  TEST_CASE("TetOverlapVolumeAABB") {
    const auto sphere = Sphere{Vector{2, 0, 0}, 0.5};

    create_benchmark("tet_overlap_volume[AABB]", [&]() {
      const auto result = overlap_volume(sphere, tet);
      ankerl::nanobench::doNotOptimizeAway(result);
    });
  }
}
