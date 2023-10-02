/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2022 Severin Strobl <severin.strobl@fau.de>
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

#include <doctest/doctest.h>
#include <nanobench.h>

#include "overlap/overlap.hpp"

TEST_CASE("NormalNewell") {
  using namespace overlap::detail;

  const auto points = std::array<Vector, 3>{
      {{-0.8482081444352685, -0.106496132943784, -0.5188463331100054},
       {-0.8482081363047198, -0.1064961977010221, -0.5188463331100054},
       {-0.8482081363047198, -0.106496132943784, -0.5188463464017972}}};

  const auto center =
      (Scalar{1} / Scalar{3}) *
      std::accumulate(points.begin(), points.end(), Vector::Zero().eval());

  auto log = std::ofstream{"normal_newell.json"};
  ankerl::nanobench::Bench()
      .title("normal_newell")
      .run("normal_newell",
           [&]() {
             const auto normal =
                 normal_newell(points.begin(), points.end(), center);
             ankerl::nanobench::doNotOptimizeAway(normal);
           })
      .render(ankerl::nanobench::templates::pyperf(), log);
}

TEST_CASE("RegularizedWedge") {
  using namespace overlap::detail;

  auto log = std::ofstream{"regularized_wedge.json"};
  auto rng = ankerl::nanobench::Rng{};
  ankerl::nanobench::Bench()
      .title("regularized_wedge")
      .run("regularized_wedge",
           [&]() {
             const auto result = regularized_wedge(1.0, rng.uniform01(),
                                                   rng.uniform01() * 0.5 * pi);

             ankerl::nanobench::doNotOptimizeAway(result);
           })
      .doNotOptimizeAway(rng)
      .render(ankerl::nanobench::templates::pyperf(), log);
  ;
}
