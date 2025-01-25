// Copyright (C) 2025 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

#include "overlap/overlap.hpp"

TEST_SUITE("line_sphere_intersection") {
  using namespace overlap;

  const auto unit_sphere = Sphere{};

  TEST_CASE("NoIntersection") {
    const auto base = Vector{2, 0, 0};
    const auto direction = Vector::UnitZ();

    for (const auto& intersection :
         detail::line_sphere_intersection(base, direction, unit_sphere)) {
      CHECK_FALSE(intersection.has_value());
    }
  }

  TEST_CASE("Tangential") {
    const auto base = Vector{1, 0, -1};
    const auto direction = 2 * Vector::UnitZ();

    const auto intersections =
        detail::line_sphere_intersection(base, direction, unit_sphere);

    CHECK_FALSE(intersections[1].has_value());
    REQUIRE(intersections[0].has_value());
    CHECK(*intersections[0] == 0.5);
  }

  TEST_CASE("Intersection") {
    const auto base = Vector{-2, 0, 0};
    const auto direction = Vector::UnitX();

    for (const auto radius : {0.5, 1.0, 2.0}) {
      const auto sphere = Sphere{Vector::Zero(), radius};
      const auto intersections =
          detail::line_sphere_intersection(base, direction, sphere);

      REQUIRE(intersections[0].has_value());
      CHECK(*intersections[0] == -base.x() - sphere.radius);

      REQUIRE(intersections[1].has_value());
      CHECK(*intersections[1] == *intersections[0] + 2.0 * sphere.radius);
    }
  }

  TEST_CASE("VanishingDirection") {
    const auto base = Vector::Zero();
    const auto direction = Vector::Zero();

    for (const auto& intersection :
         detail::line_sphere_intersection(base, direction, unit_sphere)) {
      CHECK_FALSE(intersection.has_value());
    }
  }
}
