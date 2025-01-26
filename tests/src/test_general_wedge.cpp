// Copyright (C) 2025 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

#include "overlap/overlap.hpp"

#include <limits>
#include <type_traits>

template<unsigned int N>
struct Dim : public std::integral_constant<unsigned int, N> {};

TEST_SUITE("general_wedge") {
  using namespace overlap;

  const auto unit_sphere = Sphere{};

  TEST_CASE_TEMPLATE("SimpleSphericalWedge", Dimensionality, Dim<2>, Dim<3>) {
    constexpr auto dim = Dimensionality::value;

    const auto p0 = detail::Plane{Vector::Zero(), Vector::UnitX()};
    const auto p1 = detail::Plane{Vector::Zero(), Vector::UnitY()};
    const auto d = Vector::Zero();

    CHECK_NOTHROW({
      const auto result = detail::general_wedge<dim>(unit_sphere, p0, p1, d);

      if constexpr (dim == 2u) {
        CHECK(
            result ==
            Approx(detail::pi).epsilon(std::numeric_limits<Scalar>::epsilon()));
      }

      if constexpr (dim == 3u) {
        CHECK(result == Approx((1.0 / 3.0) * detail::pi)
                            .epsilon(std::numeric_limits<Scalar>::epsilon()));
      }
    });
  }

  TEST_CASE_TEMPLATE("SingleRegularizedWedge", Dimensionality, Dim<2>, Dim<3>) {
    constexpr auto dim = Dimensionality::value;

    const auto d = 0.5 * Vector::UnitX();
    const auto p0 = detail::Plane{Vector::Zero(), Vector::UnitY()};
    const auto p1 = detail::Plane{d, Vector::UnitX()};

    CHECK_NOTHROW({
      const auto result = detail::general_wedge<dim>(unit_sphere, p0, p1, d);

      if constexpr (dim == 2u) {
        CHECK(result == Approx(0.5 * unit_sphere.cap_surface_area(0.5))
                            .epsilon(std::numeric_limits<Scalar>::epsilon()));
      }

      if constexpr (dim == 3u) {
        CHECK(result == Approx(0.5 * unit_sphere.cap_volume(0.5))
                            .epsilon(std::numeric_limits<Scalar>::epsilon()));
      }
    });
  }

  TEST_CASE_TEMPLATE("Tangential", Dimensionality, Dim<2>, Dim<3>) {
    constexpr auto dim = Dimensionality::value;

    const auto p0 = detail::Plane{Vector::Zero(), Vector::UnitX()};
    const auto p1 = detail::Plane{Vector::Zero(), Vector::UnitY()};
    const auto d = (Vector::UnitX() + Vector::UnitY()).normalized();

    CHECK_NOTHROW({
      const auto result = detail::general_wedge<dim>(unit_sphere, p0, p1, d);
      CHECK(result == Scalar{0});
    });
  }
}
