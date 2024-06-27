// Copyright (C) 2021-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("contains") {
  using namespace overlap;

  TEST_CASE("SpherePoint") {
    CHECK(detail::contains(Sphere{}, Vector::Zero()));
    CHECK(detail::contains(Sphere{}, Vector::Constant(0.25)));
    CHECK(!detail::contains(Sphere{}, Vector::Constant(2.0)));
  }

  TEST_CASE("PolygonPoint") {
    const auto tri = detail::Triangle{{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}}};

    CHECK(detail::contains(tri, Vector{0.25, 0.25, 0.0}));
    CHECK(!detail::contains(tri, Vector{1.0, 1.0, 0.0}));

    const auto quad =
        detail::Quadrilateral{{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}}};

    CHECK(detail::contains(quad, Vector{0.5, 0.5, 0.0}));
    CHECK(!detail::contains(quad, Vector{-1.0, 0.0, 0.0}));
  }

  TEST_CASE("ElementPoint") {
    CHECK(detail::contains(unit_hexahedron(), Vector::Zero()));
    CHECK(detail::contains(unit_hexahedron(), Vector::Constant(0.5)));
    CHECK(!detail::contains(unit_hexahedron(), Vector::Constant(2.0)));
  }

  TEST_CASE("SphereElement") {
    CHECK(detail::contains(Sphere{Vector::Zero(), 3.0}, unit_hexahedron()));
    CHECK(!detail::contains(Sphere{Vector::Zero(), 0.5}, unit_hexahedron()));
  }
}
