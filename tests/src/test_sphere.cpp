/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2021-2022 Severin Strobl <git@severin-strobl.de>
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

#include "common.hpp"

TEST_SUITE("Sphere") {
  using namespace overlap;

  constexpr auto epsilon = std::numeric_limits<Scalar>::epsilon();
  const auto s = Sphere{};

  TEST_CASE("Volume") {
    CHECK(s.volume == Approx(4.0 / 3.0 * detail::pi).epsilon(epsilon));
  }

  TEST_CASE("CapVolume") {
    CHECK(s.cap_volume(-1 * s.radius) == Scalar{0});
    CHECK(s.cap_volume(0) == Scalar{0});
    CHECK(s.cap_volume(0.5 * s.radius) ==
          Approx(0.625 * detail::pi / 3.0).epsilon(epsilon));
    CHECK(s.cap_volume(s.radius) == Approx(0.5 * s.volume).epsilon(epsilon));
    CHECK(s.cap_volume(2 * s.radius) == s.volume);
    CHECK(s.cap_volume(3 * s.radius) == s.volume);
  }

  TEST_CASE("SurfaceArea") {
    CHECK(s.surface_area() == Approx(4.0 * detail::pi).epsilon(epsilon));
  }

  TEST_CASE("CapSurfaceArea") {
    CHECK(s.cap_surface_area(-1 * s.radius) == Scalar{0});
    CHECK(s.cap_surface_area(0) == Scalar{0});
    CHECK(s.cap_surface_area(s.radius) == 0.5 * s.surface_area());
    CHECK(s.cap_surface_area(2 * s.radius) == s.surface_area());
    CHECK(s.cap_surface_area(3 * s.radius) == s.surface_area());
  }

  TEST_CASE("DiskArea") {
    CHECK(s.disk_area(-1 * s.radius) == Scalar{0});
    CHECK(s.disk_area(0) == Scalar{0});
    CHECK(s.disk_area(s.radius) == Approx(detail::pi).epsilon(epsilon));
    CHECK(s.disk_area(2 * s.radius) == Scalar{0});
    CHECK(s.disk_area(3 * s.radius) == Scalar{0});
  }

  TEST_CASE("Contains") {
    const auto s = Sphere{Vector::Zero(), 2.0};

    CHECK(contains(s, Vector{1, 1, 1}));
    CHECK(!contains(s, Vector{2, 2, 2}));
  }
}