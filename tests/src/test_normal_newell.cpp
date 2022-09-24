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

#include "common.hpp"

TEST_SUITE("NormalNewell") {
  inline auto format_msg(const overlap::Vector& normal,
                         const overlap::Vector& expected)
      ->std::string {
    std::stringstream strm;
    strm << "invalid normal generated: [" << normal.transpose()
         << "], expected: [" << expected.transpose() << "]";

    return strm.str();
  }

  TEST_CASE("Simple") {
    using namespace overlap::detail;

    const auto points =
        std::array<Vector, 3>{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}};
    const auto center =
        (Scalar{1} / Scalar{3}) *
        std::accumulate(points.begin(), points.end(), Vector::Zero().eval());

    const auto normal = normal_newell(points.begin(), points.end(), center);
    const auto expected = Vector::UnitZ();

    REQUIRE_MESSAGE(
        (normal - expected).norm() < std::numeric_limits<Scalar>::epsilon(),
        format_msg(normal, expected));
  }

  TEST_CASE("EdgeCase") {
    using namespace overlap::detail;

    const std::array<Vector, 3> points = {
        {{-0.8482081444352685, -0.106496132943784, -0.5188463331100054},
         {-0.8482081363047198, -0.1064961977010221, -0.5188463331100054},
         {-0.8482081363047198, -0.106496132943784, -0.5188463464017972}}};

    const auto center =
        (Scalar{1} / Scalar{3}) *
        std::accumulate(points.begin(), points.end(), Vector::Zero().eval());

    const auto normal = normal_newell(points.begin(), points.end(), center);
    const auto expected =
        Vector{0.8482081353353663, 0.1064961653160474, 0.5188463413419023};

    REQUIRE_MESSAGE(
        (normal - expected).norm() < std::numeric_limits<Scalar>::epsilon(),
        format_msg(normal, expected));
  }

  TEST_CASE("Degenerated") {
    using namespace overlap::detail;

    const std::array<Vector, 3> points = {{{0, 0, 0}, {1, 1, 0}, {0, 0, 0}}};

    const auto center =
        (Scalar{1} / Scalar{points.size()}) *
        std::accumulate(points.begin(), points.end(), Vector::Zero().eval());

    const auto normal = normal_newell(points.begin(), points.end(), center);

    REQUIRE_MESSAGE(normal.norm() < std::numeric_limits<Scalar>::epsilon(),
                    format_msg(normal, Vector::Zero()));
  }
}