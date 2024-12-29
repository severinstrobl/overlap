// Copyright (C) 2021-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("Polygon") {
  using namespace overlap::detail;

  TEST_CASE("IsPlanarTri") {
    const auto tri =
        Triangle{{{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}}};

    CHECK(tri.is_planar());
  }

  TEST_CASE("IsPlanarQuad") {
    const auto quad0 = Quadrilateral{
        {{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}}}};

    CHECK(quad0.is_planar());

    const auto quad1 = Quadrilateral{
        {{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.1}}}};

    CHECK(!quad1.is_planar());
  }
}
