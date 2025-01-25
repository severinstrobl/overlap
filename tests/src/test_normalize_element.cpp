// Copyright (C) 2025 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

#include "overlap/overlap.hpp"

TEST_SUITE("normalize_element") {
  using namespace overlap;

  TEST_CASE("Translate") {
    const auto sphere = Sphere{Vector{3, 2, 1}, 1.0};
    const auto reference = unit_hexahedron();
    const auto translated = detail::normalize_element(sphere, reference);

    REQUIRE_EQ(translated.volume, reference.volume);
    REQUIRE_EQ(translated.surface_area(), reference.surface_area());

    for (auto vertex_idx = 0u; vertex_idx < detail::num_vertices<Hexahedron>();
         ++vertex_idx) {
      REQUIRE_EQ(translated.vertices[vertex_idx],
                 reference.vertices[vertex_idx] - sphere.center);
    }

    for (auto face_idx = 0u; face_idx < detail::num_faces<Hexahedron>();
         ++face_idx) {
      REQUIRE_EQ(translated.faces[face_idx].area,
                 reference.faces[face_idx].area);
      REQUIRE_EQ(translated.faces[face_idx].center,
                 reference.faces[face_idx].center - sphere.center);
    }
  }

  TEST_CASE("Scale") {
    const auto sphere = Sphere{Vector::Zero(), 2.0};
    const auto reference = unit_hexahedron();
    const auto translated = detail::normalize_element(sphere, reference);

    REQUIRE_EQ(translated.volume, 0.125 * reference.volume);
    REQUIRE_EQ(translated.surface_area(), 0.25 * reference.surface_area());

    for (auto vertex_idx = 0u; vertex_idx < detail::num_vertices<Hexahedron>();
         ++vertex_idx) {
      REQUIRE_EQ(translated.vertices[vertex_idx],
                 0.5 * reference.vertices[vertex_idx]);
    }

    for (auto face_idx = 0u; face_idx < detail::num_faces<Hexahedron>();
         ++face_idx) {
      REQUIRE_EQ(translated.faces[face_idx].area,
                 0.25 * reference.faces[face_idx].area);
      REQUIRE_EQ(translated.faces[face_idx].center,
                 0.5 * reference.faces[face_idx].center);
    }
  }
}
