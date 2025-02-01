// Copyright (C) 2025 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

#include "overlap/overlap.hpp"

TEST_SUITE("detect_non_planar_faces") {
  using namespace overlap;

  TEST_CASE("Planar") {
    CHECK_NOTHROW(detail::detect_non_planar_faces(unit_hexahedron()));
  }

  TEST_CASE("NonPlanar") {
    auto vertices = unit_hexahedron().vertices;
    for (auto vertex_idx = 0u; vertex_idx < 7u; ++vertex_idx) {
      vertices[vertex_idx] += Vector{0.001, 0.001, 0.001};
      if (vertex_idx == 3u) {  // valid element with enlarged base
        continue;
      }

      CHECK_THROWS_WITH_AS(
          detail::detect_non_planar_faces(Hexahedron{vertices}),
          doctest::Contains("non-planer face detected in element"),
          std::invalid_argument);
    }
  }
}
