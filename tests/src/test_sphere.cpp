// Copyright (C) 2021-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

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
