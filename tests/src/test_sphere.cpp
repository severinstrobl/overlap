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
  const auto sphere = Sphere{};

  TEST_CASE("Volume") {
    CHECK(sphere.volume == Approx(4.0 / 3.0 * detail::pi).epsilon(epsilon));
  }

  TEST_CASE("CapVolume") {
    CHECK(sphere.cap_volume(-1 * sphere.radius) == Scalar{0});
    CHECK(sphere.cap_volume(0) == Scalar{0});
    CHECK(sphere.cap_volume(0.5 * sphere.radius) ==
          Approx(0.625 * detail::pi / 3.0).epsilon(epsilon));

    CHECK(sphere.cap_volume(sphere.radius) ==
          Approx(0.5 * sphere.volume).epsilon(epsilon));

    CHECK(sphere.cap_volume(2 * sphere.radius) == sphere.volume);
    CHECK(sphere.cap_volume(3 * sphere.radius) == sphere.volume);
  }

  TEST_CASE("SurfaceArea") {
    CHECK(sphere.surface_area() == Approx(4.0 * detail::pi).epsilon(epsilon));
  }

  TEST_CASE("CapSurfaceArea") {
    CHECK(sphere.cap_surface_area(-1 * sphere.radius) == Scalar{0});
    CHECK(sphere.cap_surface_area(0) == Scalar{0});
    CHECK(sphere.cap_surface_area(sphere.radius) ==
          0.5 * sphere.surface_area());

    CHECK(sphere.cap_surface_area(2 * sphere.radius) == sphere.surface_area());
    CHECK(sphere.cap_surface_area(3 * sphere.radius) == sphere.surface_area());
  }

  TEST_CASE("DiskArea") {
    CHECK(sphere.disk_area(-1 * sphere.radius) == Scalar{0});
    CHECK(sphere.disk_area(0) == Scalar{0});
    CHECK(sphere.disk_area(sphere.radius) ==
          Approx(detail::pi).epsilon(epsilon));

    CHECK(sphere.disk_area(2 * sphere.radius) == Scalar{0});
    CHECK(sphere.disk_area(3 * sphere.radius) == Scalar{0});
  }

  TEST_CASE("Contains") {
    const auto s = Sphere{Vector::Zero(), 2.0};

    CHECK(contains(s, Vector{1, 1, 1}));
    CHECK(!contains(s, Vector{2, 2, 2}));
  }
}