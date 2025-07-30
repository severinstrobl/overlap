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
    const auto sphere = Sphere{Vector::Zero(), 1.0};

    create_benchmark("hex_overlap_volume[sphere-in-hex]", [&]() {
      const auto result = overlap_volume(sphere, hex);
      ankerl::nanobench::doNotOptimizeAway(result);
    });
  }

  TEST_CASE("HexOverlapVolume") {
    const auto sphere = Sphere{Vector::Zero(), 5.0};

    create_benchmark("hex_overlap_volume[hex-in-sphere]", [&]() {
      const auto result = overlap_volume(sphere, hex);
      ankerl::nanobench::doNotOptimizeAway(result);
    });
  }

  TEST_CASE("HexOverlapVolumeAABB") {
    const auto sphere = Sphere{Vector{5, 0, 0}, 1};

    create_benchmark("hex_overlap_volume[AABB]", [&]() {
      const auto result = overlap_volume(sphere, hex);
      ankerl::nanobench::doNotOptimizeAway(result);
    });
  }

  TEST_CASE("HexOverlapVolumeRandomized") {
    constexpr auto seed = 79'866'982'766'580U;
    auto rng = ankerl::nanobench::Rng{seed};

    create_benchmark("hex_overlap_volume[random]", [&rng]() {
      const auto radius = (2.4 * rng.uniform01()) + 0.1;
      const auto center =
          4.0 * Vector{rng.uniform01(), rng.uniform01(), rng.uniform01()} -
          Vector::Constant(2.0);

      const auto sphere = overlap::Sphere{center, radius};
      const auto result = overlap_volume(sphere, hex);
      ankerl::nanobench::doNotOptimizeAway(result);
    }).epochIterations(25'000);
  }
}
