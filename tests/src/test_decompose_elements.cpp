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

#include "gtest/gtest.h"

#include "overlap/overlap.hpp"

#include "common.hpp"

using namespace overlap;

TEST(DecomposeElements, Hexahedron) {
  Hexahedron hex = unit_hexahedron();

  std::array<Tetrahedron, 4> subTets;
  std::array<Tetrahedron, 5> tets5;
  std::array<Tetrahedron, 6> tets6;
  std::array<Wedge, 2> wedges;

  decompose(hex, tets5);
  decompose(hex, tets6);
  decompose(hex, wedges);

  auto tets5Volume = Scalar{0};
  for (const auto& tet : tets5) {
    tets5Volume += tet.volume;
  }

  auto tets6Volume = Scalar{0};
  auto tets24Volume = Scalar{0};
  for (const auto& tet : tets6) {
    decompose(tet, subTets);
    tets6Volume += tet.volume;

    for (const auto& subTet : subTets) {
      tets24Volume += subTet.volume;
    }
  }

  constexpr auto delta = Scalar{5e2} * std::numeric_limits<Scalar>::epsilon();

  ASSERT_NEAR(hex.volume, tets5Volume, delta);
  ASSERT_NEAR(hex.volume, tets6Volume, delta);
  ASSERT_NEAR(hex.volume, tets24Volume, delta);
  ASSERT_NEAR(hex.volume, wedges[0].volume + wedges[1].volume, delta);
}
