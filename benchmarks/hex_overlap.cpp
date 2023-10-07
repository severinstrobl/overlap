/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2022 Severin Strobl <severin.strobl@fau.de>
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
