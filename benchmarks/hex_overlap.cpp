// Copyright (C) 2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include <doctest/doctest.h>

#include "overlap/overlap.hpp"

#include "common.hpp"

TEST_SUITE("HexOverlap") {
  using namespace overlap;

  // clang-format off
  const auto hex = Hexahedron{{{
      {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
      {-1, -1,  1}, {1, -1,  1}, {1, 1,  1}, {-1, 1,  1},
  }}};
  // clang-format on

  TEST_CASE("HexOverlapVolume") {
    const auto sphere = Sphere{Vector::Zero(), 1.5};

    create_benchmark("hex_overlap_volume", [&]() {
      const auto result = overlap_volume(sphere, hex);
      ankerl::nanobench::doNotOptimizeAway(result);
    });
  }

  TEST_CASE("HexOverlapVolumeAABB") {
    const auto sphere = Sphere{Vector{5, 0, 0}, 1};

    create_benchmark("hex_overlap_volume_aabb", [&]() {
      const auto result = overlap_volume(sphere, hex);
      ankerl::nanobench::doNotOptimizeAway(result);
    });
  }
}
