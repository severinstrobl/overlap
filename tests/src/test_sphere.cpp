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

using namespace overlap;

TEST(Sphere, Volume) {
  const auto s = Sphere{};

  constexpr auto eps = std::numeric_limits<Scalar>::epsilon();
  ASSERT_NEAR(s.volume, 4.0 / 3.0 * detail::pi, eps);

  ASSERT_EQ(s.cap_volume(-1 * s.radius), 0.0);
  ASSERT_EQ(s.cap_volume(0), 0.0);
  ASSERT_NEAR(s.cap_volume(0.5 * s.radius), (detail::pi * 0.25 / 3.0) * 2.5,
              eps);
  ASSERT_NEAR(s.cap_volume(s.radius), 0.5 * s.volume, eps);
  ASSERT_EQ(s.cap_volume(2 * s.radius), s.volume);
  ASSERT_EQ(s.cap_volume(3 * s.radius), s.volume);
}

TEST(Sphere, SurfaceArea) {
  const auto s = Sphere{};

  constexpr auto eps = std::numeric_limits<Scalar>::epsilon();

  ASSERT_NEAR(s.surface_area(), 4.0 * detail::pi, eps);

  ASSERT_EQ(s.cap_surface_area(-1 * s.radius), 0.0);
  ASSERT_EQ(s.cap_surface_area(0), 0.0);
  ASSERT_EQ(s.cap_surface_area(s.radius), 0.5 * s.surface_area());
  ASSERT_EQ(s.cap_surface_area(2 * s.radius), s.surface_area());
  ASSERT_EQ(s.cap_surface_area(3 * s.radius), s.surface_area());
}

TEST(Sphere, DiskArea) {
  const auto s = Sphere{};

  constexpr auto eps = std::numeric_limits<Scalar>::epsilon();

  ASSERT_EQ(s.disk_area(-1 * s.radius), 0.0);
  ASSERT_EQ(s.disk_area(0), 0.0);
  ASSERT_NEAR(s.disk_area(s.radius), detail::pi, eps);
  ASSERT_EQ(s.disk_area(2 * s.radius), 0.0);
  ASSERT_EQ(s.disk_area(3 * s.radius), 0.0);
}

TEST(Sphere, Contains) {
  const auto s = Sphere{Vector::Zero(), 2.0};

  ASSERT_TRUE(contains(s, Vector{1, 1, 1}));
  ASSERT_FALSE(contains(s, Vector{2, 2, 2}));
}
