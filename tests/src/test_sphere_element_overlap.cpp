// Copyright (C) 2015-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("SphereElementOverlap") {
  using namespace overlap;

  static const auto epsilon = std::sqrt(std::numeric_limits<Scalar>::epsilon());

  // Sphere outside of hexahedron, touching one face.
  TEST_CASE("Face") {
    const auto sphere = Sphere{{0, 2, 0}, 1};

    validate_overlap_volume(sphere, unit_hexahedron(), epsilon, Scalar{0});
  }

  // Sphere outside of hexahedron, intersecting one face, touching 4 edges
  TEST_CASE("FaceMaxOverlap") {
    const auto sphere = Sphere{{1, 0, 0}, 1};

    validate_overlap_volume(sphere, unit_hexahedron(), epsilon,
                            0.5 * sphere.volume);
  }

  // Sphere intersects one edge (and thus 2 faces).
  TEST_CASE("Edge") {
    const auto sphere = Sphere{Vector{0, -1, 1}, 1};

    validate_overlap_volume(sphere, unit_hexahedron(), epsilon,
                            0.25 * sphere.volume);
  }

  // Sphere outside of hexahedron, touching one edge.
  TEST_CASE("EdgeTouching") {
    const auto offset =
        1.0 + 0.5 * std::sqrt(2.0) - Scalar{1e2} * detail::tiny_epsilon;
    const auto sphere = Sphere{Vector{offset, offset, 0}, 1};

    validate_overlap_volume(sphere, unit_hexahedron(), epsilon, Scalar{0});
  }

  // Sphere outside of hexahedron, touching one edge, with the sphere centered
  // above the edge and shifted along the edge.
  TEST_CASE("EdgeTouchingCentered") {
    const auto radius = 0.05;
    const auto offsets = {0.0, 1e-6, 0.005, 0.01};

    for (const auto offset : offsets) {
      CAPTURE(offset);

      const auto sphere = Sphere{Vector{1.0, 1.0 + radius, offset}, radius};
      validate_overlap_volume(sphere, unit_hexahedron(), epsilon, Scalar{0});
    }
  }

  // Sphere intersects one edge (and thus 2 faces), edge passing through
  // center of sphere -> spherical wedge with angle pi/4.
  TEST_CASE("Wedge") {
    const auto sphere = Sphere{Vector::Zero(), 1};

    // clang-format off
    const auto hex = Hexahedron{{{
      {0, 0, -1}, {2, 2, -1},
      {2, 4, -1}, {0, 4, -1},

      {0, 0, 1}, {2, 2, 1},
      {2, 4, 1}, {0, 4, 1}}}};
    // clang-format on

    validate_overlap_volume(
        sphere, hex, epsilon,
        2.0 / 3.0 * sphere.radius * sphere.radius * 0.25 * detail::pi);
  }

  // Sphere intersects one vertex (and thus 3 edges and 3 faces)
  TEST_CASE("Vertex") {
    const auto sphere = Sphere{Vector{1, -1, 1}, 1};

    validate_overlap_volume(sphere, unit_hexahedron(), epsilon,
                            0.125 * sphere.volume);
  }

  // Sphere outside of hexahedron, touching one vertex.
  TEST_CASE("VertexTouching") {
    const auto sphere = Sphere{Vector{2, -1, 1}, 1};

    validate_overlap_volume(sphere, unit_hexahedron(), epsilon, Scalar{0});
  }

  // Sphere outside of hexahedron, slightly overlapping one vertex.
  TEST_CASE("VertexOverlap") {
    const auto sphere = Sphere{Vector{2 - 10 * detail::tiny_epsilon, -1, 1}, 1};

    validate_overlap_volume(sphere, unit_hexahedron(), epsilon, Scalar{0});
  }

  // Sphere outside of hexahedron, slightly overlapping one vertex.
  TEST_CASE("VertexTinyOverlap") {
    const auto offset =
        1.0 + 0.5 * std::sqrt(2.0) - Scalar{1e2} * detail::tiny_epsilon;
    const auto sphere = Sphere{Vector{offset, offset, 1}, 1};

    validate_overlap_volume(sphere, unit_hexahedron(), epsilon, Scalar{0});
  }

  // Sphere contains hexahedron.
  TEST_CASE("HexInSphere") {
    const auto sphere = Sphere{Vector::Zero(), 2.0};
    const auto hex = unit_hexahedron();

    validate_overlap_volume(sphere, hex, epsilon, hex.volume);
  }

  // Sphere contained in hexahedron.
  TEST_CASE("SphereInHex") {
    const auto sphere = Sphere{Vector::Zero(), 0.5};

    validate_overlap_volume(sphere, unit_hexahedron(), epsilon, sphere.volume);
  }

  // Ensure non-planar faces are detected.
  TEST_CASE("NonPlanarFaces") {
    auto vertices = unit_hexahedron().vertices;
    vertices[0] += Vector{0, 0, -0.25};

    REQUIRE_THROWS_AS(overlap_volume(Sphere{}, Hexahedron{vertices}),
                      std::invalid_argument);
  }
}
