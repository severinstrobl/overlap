// Copyright (C) 2025 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

#include "overlap/overlap.hpp"

#include <utility>
#include <vector>

TEST_SUITE("unit_sphere_intersections") {
  using namespace overlap;

  TEST_CASE("NoIntersection") {
    const auto&& [entity_intersections, edge_intersections] =
        unit_sphere_intersections(unit_hexahedron());

    CHECK_FALSE(entity_intersections.vertices.any());
    CHECK_FALSE(entity_intersections.edges.any());
    CHECK_FALSE(entity_intersections.faces.any());

    std::ignore = edge_intersections;
  }

  TEST_CASE("FaceIntersection") {
    const auto shifts = std::vector<Vector>{Vector::UnitZ(),  Vector::UnitY(),
                                            -Vector::UnitX(), -Vector::UnitY(),
                                            Vector::UnitX(),  -Vector::UnitZ()};

    for (auto face_idx = 0u; face_idx < detail::num_faces<Hexahedron>();
         ++face_idx) {
      const auto transformation =
          detail::Transformation{shifts[face_idx], Scalar{1}};

      auto hex = unit_hexahedron();
      hex.apply(transformation);

      const auto&& [entity_intersections, edge_intersections] =
          unit_sphere_intersections(hex);

      CHECK_FALSE(entity_intersections.vertices.any());
      CHECK_FALSE(entity_intersections.edges.any());

      CHECK_EQ(entity_intersections.faces.count(), 1u);
      CHECK(entity_intersections.faces[face_idx]);

      std::ignore = edge_intersections;
    }
  }

  // TODO: edge and vertex intersections
}
