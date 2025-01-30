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

    for (const auto& intersection : edge_intersections) {
      CHECK_FALSE(intersection.has_value());
    }
  }

  TEST_CASE("FaceIntersection") {
    const auto shifts = std::array<Vector, detail::num_faces<Hexahedron>()>{
        Vector::UnitZ(),  Vector::UnitY(), -Vector::UnitX(),
        -Vector::UnitY(), Vector::UnitX(), -Vector::UnitZ(),
    };

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

      for (const auto& intersection : edge_intersections) {
        CHECK_FALSE(intersection.has_value());
      }
    }
  }

  TEST_CASE("EdgeIntersection") {
    const auto shifts = std::array<Vector, detail::num_edges<Hexahedron>()>{
        Vector::UnitY() + Vector::UnitZ(),  -Vector::UnitX() + Vector::UnitZ(),
        -Vector::UnitY() + Vector::UnitZ(), Vector::UnitX() + Vector::UnitZ(),
        Vector::UnitX() + Vector::UnitY(),  -Vector::UnitX() + Vector::UnitY(),
        -Vector::UnitX() - Vector::UnitY(), Vector::UnitX() - Vector::UnitY(),
        Vector::UnitY() - Vector::UnitZ(),  -Vector::UnitX() - Vector::UnitZ(),
        -Vector::UnitY() - Vector::UnitZ(), Vector::UnitX() - Vector::UnitZ(),
    };

    const auto faces = std::array<std::array<std::uint8_t, 2>,
                                  detail::num_edges<Hexahedron>()>{{
        {{0u, 1u}},
        {{0u, 2u}},
        {{0u, 3u}},
        {{0u, 4u}},
        {{1u, 4u}},
        {{1u, 2u}},
        {{2u, 3u}},
        {{3u, 4u}},
        {{1u, 5u}},
        {{2u, 5u}},
        {{3u, 5u}},
        {{4u, 5u}},
    }};

    const auto cwiseNegation = [](const Vector& v) {
      return (Vector::Ones() - v.cwiseAbs()).eval();
    };

    for (auto edge_idx = 0u; edge_idx < detail::num_edges<Hexahedron>();
         ++edge_idx) {
      const auto transformation =
          detail::Transformation{shifts[edge_idx], Scalar{2}};

      auto hex = unit_hexahedron();
      hex.apply(transformation);

      const auto&& [entity_intersections, edge_intersections] =
          unit_sphere_intersections(hex);

      CHECK_FALSE(entity_intersections.vertices.any());

      CHECK_EQ(entity_intersections.edges.count(), 1u);
      CHECK(entity_intersections.edges[edge_idx]);

      CHECK_EQ(entity_intersections.faces.count(), 2u);
      for (const auto face_idx : faces[edge_idx]) {
        CHECK(entity_intersections.faces[face_idx]);
      }

      for (auto intersect_edge_idx = 0u;
           intersect_edge_idx < detail::num_edges<Hexahedron>();
           ++intersect_edge_idx) {
        if (intersect_edge_idx == edge_idx) {
          REQUIRE(edge_intersections[intersect_edge_idx].has_value());

          const auto inverse =
              (edge_intersections[intersect_edge_idx].value()[0].array() <
               Scalar{0})
                  .any();

          const auto offset = (cwiseNegation(shifts[edge_idx]) *
                               (inverse ? Scalar{-1} : Scalar{1}))
                                  .eval();

          CHECK_FALSE(
              (edge_intersections[intersect_edge_idx].value()[0].array() -
               offset.array())
                  .abs()
                  .any());

          CHECK_FALSE(
              (edge_intersections[intersect_edge_idx].value()[1].array() +
               offset.array())
                  .abs()
                  .any());
        } else {
          CHECK_FALSE(edge_intersections[intersect_edge_idx].has_value());
        }
      }
    }
  }

  TEST_CASE("VertexIntersection") {
    const auto shifts = std::array<Vector, detail::num_edges<Hexahedron>()>{
        Vector::UnitX() + Vector::UnitY() + Vector::UnitZ(),
        -Vector::UnitX() + Vector::UnitY() + Vector::UnitZ(),
        -Vector::UnitX() - Vector::UnitY() + Vector::UnitZ(),
        Vector::UnitX() - Vector::UnitY() + Vector::UnitZ(),
        Vector::UnitX() + Vector::UnitY() - Vector::UnitZ(),
        -Vector::UnitX() + Vector::UnitY() - Vector::UnitZ(),
        -Vector::UnitX() - Vector::UnitY() - Vector::UnitZ(),
        Vector::UnitX() - Vector::UnitY() - Vector::UnitZ(),
    };

    const auto faces = std::array<std::array<std::uint8_t, 3>,
                                  detail::num_vertices<Hexahedron>()>{{
        {{0u, 1u, 4u}},
        {{0u, 1u, 2u}},
        {{0u, 2u, 3u}},
        {{0u, 3u, 4u}},
        {{1u, 4u, 5u}},
        {{1u, 2u, 5u}},
        {{2u, 3u, 5u}},
        {{3u, 4u, 5u}},
    }};

    const auto edges = std::array<std::array<std::uint8_t, 3>,
                                  detail::num_vertices<Hexahedron>()>{{
        {{0u, 3u, 4u}},
        {{0u, 1u, 5u}},
        {{2u, 2u, 6u}},
        {{2u, 3u, 7u}},
        {{4u, 8u, 11u}},
        {{5u, 8u, 9u}},
        {{6u, 9u, 10u}},
        {{7u, 10u, 11u}},
    }};

    for (auto vertex_idx = 0u; vertex_idx < detail::num_vertices<Hexahedron>();
         ++vertex_idx) {
      const auto transformation =
          detail::Transformation{shifts[vertex_idx], Scalar{1}};

      auto hex = unit_hexahedron();
      hex.apply(transformation);

      const auto&& [entity_intersections, edge_intersections] =
          unit_sphere_intersections(hex);

      CHECK_EQ(entity_intersections.vertices.count(), 1u);
      CHECK(entity_intersections.vertices[vertex_idx]);

      CHECK_EQ(entity_intersections.edges.count(), 3u);
      for (const auto edge_idx : edges[vertex_idx]) {
        CHECK(entity_intersections.edges[edge_idx]);

        CHECK(edge_intersections[edge_idx].has_value());
      }

      CHECK_EQ(entity_intersections.faces.count(), 3u);
      for (const auto face_idx : faces[vertex_idx]) {
        CHECK(entity_intersections.faces[face_idx]);
      }
    }
  }
}
