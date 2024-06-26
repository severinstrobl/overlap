// Copyright (C) 2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include <doctest/doctest.h>

#include "overlap/overlap.hpp"

#include "common.hpp"

TEST_CASE("NormalNewell") {
  using namespace overlap::detail;

  const auto points = std::array<Vector, 3>{
      {{-0.8482081444352685, -0.106496132943784, -0.5188463331100054},
       {-0.8482081363047198, -0.1064961977010221, -0.5188463331100054},
       {-0.8482081363047198, -0.106496132943784, -0.5188463464017972}}};

  const auto center =
      (Scalar{1} / Scalar{3}) *
      std::accumulate(points.begin(), points.end(), Vector::Zero().eval());

  create_benchmark("normal_newell", [&]() {
    const auto normal = normal_newell(points.begin(), points.end(), center);
    ankerl::nanobench::doNotOptimizeAway(normal);
  });
}

TEST_CASE("RegularizedWedge") {
  using namespace overlap::detail;

  auto rng = ankerl::nanobench::Rng{};
  create_benchmark("regularized_wedge", [&]() {
    const auto result =
        regularized_wedge(1.0, rng.uniform01(), rng.uniform01() * 0.5 * pi);
    ankerl::nanobench::doNotOptimizeAway(result);
  }).doNotOptimizeAway(rng);
}
