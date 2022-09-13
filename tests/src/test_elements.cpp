/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2021-2022 Severin Strobl <git@severin-strobl.de>
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

TEST(Tetrahedron, Constructors) {
  using namespace overlap;

  const auto tet0 = Tetrahedron{};
  ASSERT_EQ(tet0.volume, Scalar{0});

  const auto sqrt3 = std::sqrt(Scalar{3});
  const auto sqrt6 = std::sqrt(Scalar{6});

  // clang-format off
  const auto tet1 = Tetrahedron{{{
    {-sqrt3 / 6.0, -1.0 / 2.0, 0}, {sqrt3 / 3.0, 0, 0},
    {-sqrt3 / 6.0, +1.0 / 2.0, 0}, {0, 0, sqrt6 / 3.0}}}};

  const auto tet2 = Tetrahedron{std::array<Vector, 4>{{
    {-sqrt3 / 6.0, -1.0 / 2.0, 0}, {sqrt3 / 3.0, 0, 0},
    {-sqrt3 / 6.0, +1.0 / 2.0, 0}, {0, 0, sqrt6 / 3.0}}}};
  // clang-format on

  constexpr auto eps = std::numeric_limits<Scalar>::epsilon();
  ASSERT_NEAR(tet1.volume, Scalar{(1.0 / 12.0) * std::sqrt(2.0)}, eps);
  ASSERT_NEAR(tet1.volume, tet2.volume, eps);
}

TEST(Wedge, Constructors) {
  using namespace overlap;

  const auto wedge0 = Wedge{};
  ASSERT_EQ(wedge0.volume, 0.0);

  // clang-format off
  const auto wedge1 = Wedge{{{
    {-1, -1, -1}, {1, -1, -1}, {1, 1, -1},
    {-1, -1, 1},  {1, -1, 1},  {1, 1, 1}}}};

  const auto wedge2 = Wedge{std::array<Vector, 6>{{
    {-1, -1, -1}, {1, -1, -1}, {1, 1, -1},
    {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}}}};
  // clang-format on

  constexpr auto eps = std::numeric_limits<Scalar>::epsilon();
  ASSERT_NEAR(wedge1.volume, Scalar{4}, eps);
  ASSERT_NEAR(wedge1.volume, wedge2.volume, eps);
}

TEST(Hexahedron, Constructors) {
  using namespace overlap;

  const auto hex0 = Hexahedron{};
  ASSERT_EQ(hex0.volume, 0.0);

  // clang-format off
  const auto hex1 = Hexahedron{{{
    {-1, -1, -1}, {1, -1, -1},
    { 1,  1, -1}, {-1,  1, -1},

    {-1, -1,  1}, {1, -1,  1},
    { 1,  1,  1}, {-1,  1,  1}}}};

  const auto hex2 = Hexahedron{std::array<Vector, 8>{{
    {-1, -1, -1}, {1, -1, -1}, {1,  1, -1}, {-1,  1, -1},
    {-1, -1,  1}, {1, -1,  1}, {1,  1,  1}, {-1,  1,  1}}}};
  // clang-format on

  constexpr auto eps = std::numeric_limits<Scalar>::epsilon();
  ASSERT_NEAR(hex1.volume, Scalar{8}, eps);
  ASSERT_NEAR(hex1.volume, hex2.volume, eps);
}
