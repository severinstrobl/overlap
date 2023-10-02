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
#include <nanobench.h>

#include "overlap/overlap.hpp"

TEST_CASE("HexOverlapVolume") {
  using namespace overlap;

  const auto sphere = Sphere{Vector::Zero(), 1};

  // clang-format off
  const auto hex = Hexahedron{{{
    {0, 0, -1}, {2, 2, -1},
    {2, 4, -1}, {0, 4, -1},

    {0, 0, 1}, {2, 2, 1},
    {2, 4, 1}, {0, 4, 1}}}};
  // clang-format on

  auto log = std::ofstream{"hex_overlap_volume.json"};
  ankerl::nanobench::Bench()
      .title("hex_overlap_volume")
      .run("hex_overlap_volume",
           [&]() {
             const auto result = overlap_volume(sphere, hex);
             ankerl::nanobench::doNotOptimizeAway(result);
           })
      .render(ankerl::nanobench::templates::pyperf(), log);
}

TEST_CASE("HexOverlapVolumeAABB") {
  using namespace overlap;

  const auto sphere = Sphere{Vector::Zero(), 1};

  // clang-format off
  const auto hex = Hexahedron{{{
    {2, 0, 0}, {4, 0, 0},
    {4, 2, 0}, {2, 2, 0},

    {2, 0, 2}, {4, 0, 2},
    {4, 2, 2}, {2, 2, 2}}}};
  // clang-format on

  auto log = std::ofstream{"hex_overlap_volume_aabb.json"};
  ankerl::nanobench::Bench()
      .title("hex_overlap_volume_aabb")
      .run("hex_overlap_volume_aabb",
           [&]() {
             const auto result = overlap_volume(sphere, hex);
             ankerl::nanobench::doNotOptimizeAway(result);
           })
      .render(ankerl::nanobench::templates::pyperf(), log);
}
