// Copyright (C) 2015-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("NormalNewell") {
  using Scalar = overlap::Scalar;
  using Vector = overlap::Vector;

  inline auto format_msg(const Vector& normal, const Vector& expected)
      -> std::string {
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

    REQUIRE_MESSAGE((normal.norm() < std::numeric_limits<Scalar>::epsilon() ||
                     normal == Vector::UnitZ()),
                    format_msg(normal, Vector::Zero()).append([]() {
                      std::stringstream strm;
                      strm << " or [" << Vector::UnitZ().transpose() << "]";
                      return strm.str();
                    }()));
  }
}
