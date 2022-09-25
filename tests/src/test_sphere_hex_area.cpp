/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2015-2022 Severin Strobl <severin.strobl@fau.de>
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

TEST_SUITE("SphereHexAreaTest") {
  using namespace overlap;

  // Sphere intersects one face.
  TEST_CASE("Face") {
    const auto hex = unit_hexahedron();
    const auto s = Sphere{{0, 0, 1}, 0.75};

    auto result = overlap_area(s, hex);

    std::array<Scalar, 8> resultExact{};
    resultExact[0] = Scalar{0.5} * s.surface_area();
    resultExact[6] = s.disk_area(s.radius);
    resultExact[7] = resultExact[6];

    for (std::size_t i = 0; i < resultExact.size(); ++i) {
      CHECK(result[i] == Approx(resultExact[i])
                             .epsilon(std::numeric_limits<Scalar>::epsilon()));
    }
  }

  // Sphere intersects one edge (and thus 1 edge and 2 faces).
  TEST_CASE("Edge") {
    const auto hex = unit_hexahedron();
    const auto sphere = Sphere{{1, 1, 0}, 0.75};

    auto result = overlap_area(sphere, hex);

    std::array<Scalar, 8> resultExact{};
    resultExact[0] = Scalar{0.25} * sphere.surface_area();
    resultExact[3] = 0.5 * sphere.disk_area(sphere.radius);
    resultExact[4] = resultExact[3];
    resultExact[7] = 2 * resultExact[3];

    for (std::size_t i = 0; i < resultExact.size(); ++i) {
      CHECK(result[i] == Approx(resultExact[i])
                             .epsilon(std::numeric_limits<Scalar>::epsilon()));
    }
  }

  // Sphere intersects one vertex (and thus 3 edge and 3 faces).
  TEST_CASE("Vertex") {
    const auto hex = unit_hexahedron();
    Sphere s({1, 1, 1}, 0.75);

    auto result = overlap_area(s, hex);

    std::array<Scalar, 8> resultExact{};
    resultExact[0] = Scalar{0.125} * s.surface_area();
    resultExact[3] = 0.25 * s.disk_area(s.radius);
    resultExact[4] = resultExact[3];
    resultExact[6] = resultExact[3];
    resultExact[7] = 3 * resultExact[3];

    constexpr auto epsilon =
        Scalar{1e3} * std::numeric_limits<Scalar>::epsilon();

    for (std::size_t i = 0; i < resultExact.size(); ++i) {
      CHECK(result[i] == Approx(resultExact[i]).epsilon(epsilon));
    }
  }
}