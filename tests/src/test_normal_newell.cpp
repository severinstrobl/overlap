/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2015-2022 Severin Strobl <git@severin-strobl.de>
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

#include <array>
#include "common.hpp"
#include "overlap/overlap.hpp"

TEST_SUITE("NormalNewell") {
  using Scalar = overlap::Scalar;
  using Vector = overlap::Vector;

  inline auto format_msg(const Vector& normal, const Vector& expected)
      ->std::string {
    std::stringstream strm;
    strm << "invalid normal generated: [" << normal.transpose()
         << "], expected: [" << expected.transpose() << "]";

    return strm.str();
  }

  template<typename Iterator>
  inline auto calc_center(Iterator first, Iterator last) {
    return ((Scalar{1} / static_cast<Scalar>(std::distance(first, last))) *
            std::accumulate(first, last, Vector::Zero().eval()))
        .eval();
  }

  TEST_CASE("Simple") {
    const auto points =
        std::array<Vector, 3>{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}};

    const auto center = calc_center(points.begin(), points.end());
    const auto normal =
        overlap::detail::normal_newell(points.begin(), points.end(), center);

    const auto expected = Vector::UnitZ();

    REQUIRE_MESSAGE(
        (normal - expected).norm() < std::numeric_limits<Scalar>::epsilon(),
        format_msg(normal, expected));
  }

  TEST_CASE("EdgeCase") {
    const std::array<Vector, 3> points = {
        {{-0.8482081444352685, -0.106496132943784, -0.5188463331100054},
         {-0.8482081363047198, -0.1064961977010221, -0.5188463331100054},
         {-0.8482081363047198, -0.106496132943784, -0.5188463464017972}}};

    const auto center = calc_center(points.begin(), points.end());
    const auto normal =
        overlap::detail::normal_newell(points.begin(), points.end(), center);

    const auto expected =
        Vector{0.8482081353353663, 0.1064961653160474, 0.5188463413419023};

    REQUIRE_MESSAGE(
        (normal - expected).norm() < std::numeric_limits<Scalar>::epsilon(),
        format_msg(normal, expected));
  }

  TEST_CASE("Degenerated") {
    const std::array<Vector, 3> points = {{{0, 0, 0}, {1, 1, 0}, {0, 0, 0}}};

    const auto center = calc_center(points.begin(), points.end());
    const auto normal =
        overlap::detail::normal_newell(points.begin(), points.end(), center);

    REQUIRE_MESSAGE(normal.norm() < std::numeric_limits<Scalar>::epsilon(),
                    format_msg(normal, Vector::Zero()));
  }

  TEST_CASE("StarShapedPolygon") {
    for (auto i = 1U; i < 5U; ++i) {
      const auto k = static_cast<Scalar>(i);
      const auto vertices = std::array<Vector, 8>{{
          {1, 1, 1},
          {6, 6 - k, 1},
          {11, 1, 1},
          {6 + k, 6, 1},
          {11, 11, 1},
          {6, 6 + k, 1},
          {1, 11, 1},
          {6 - k, 6, 1},
      }};

      REQUIRE_EQ(overlap::detail::normal_newell(
                     std::begin(vertices), std::end(vertices),
                     calc_center(std::begin(vertices), std::end(vertices))),
                 Vector::UnitZ());
    }
  }
}
