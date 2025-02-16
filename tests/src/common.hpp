// Copyright (C) 2015-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#ifndef OVERLAP_TEST_COMMON_HPP
#define OVERLAP_TEST_COMMON_HPP

#include <iostream>
#include <string_view>

#include <doctest/doctest.h>

class AssertionError : public std::runtime_error {
 public:
  static constexpr auto assertion_failed_msg = "[overlap] assertion failed: ";

  explicit AssertionError(std::string_view msg) :
      std::runtime_error{assertion_failed_msg + std::string{msg}} {}
};

[[noreturn]] inline void throw_assertion_error(std::string_view msg) {
  throw AssertionError{msg};
}

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage,readability-identifier-naming)
#define overlap_assert(expr, msg) \
  (static_cast<bool>(expr) ? static_cast<void>(0) : throw_assertion_error(msg))

#include "overlap/overlap.hpp"

using Approx = doctest::Approx;

namespace overlap {

// constexpr version of abs, constexpr in <cmath> only with C++23
template<typename T, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
constexpr auto abs(const T value) -> T {
  return value < T{0} ? -value : value;
}

inline auto unit_hexahedron(const Scalar scaling = Scalar{1}) -> Hexahedron {
  const auto v0 = Vector{-1, -1, -1};
  const auto v1 = Vector{1, -1, -1};
  const auto v2 = Vector{1, 1, -1};
  const auto v3 = Vector{-1, 1, -1};
  const auto v4 = Vector{-1, -1, 1};
  const auto v5 = Vector{1, -1, 1};
  const auto v6 = Vector{1, 1, 1};
  const auto v7 = Vector{-1, 1, 1};

  return Hexahedron{v0 * scaling, v1 * scaling, v2 * scaling, v3 * scaling,
                    v4 * scaling, v5 * scaling, v6 * scaling, v7 * scaling};
}

inline void validate_overlap_volume(const Sphere& s, const Hexahedron& hex,
                                    const Scalar epsilon,
                                    const Scalar exactResult = Scalar{-1}) {
  CHECK_NOTHROW({
    std::array<Tetrahedron, 5> tets5;
    std::array<Tetrahedron, 6> tets6;
    std::array<Wedge, 2> wedges;

    decompose(hex, tets5);
    decompose(hex, tets6);
    decompose(hex, wedges);

    const auto overlapCalcHex = overlap_volume(s, hex);
    if (exactResult != Scalar{-1}) {
      CHECK(overlapCalcHex == Approx(exactResult).epsilon(epsilon));
    }

    const auto overlapCalcTets5 =
        overlap_volume(s, std::begin(tets5), std::end(tets5));

    const auto overlapCalcTets6 =
        overlap_volume(s, std::begin(tets6), std::end(tets6));

    auto overlapCalcTets24 = Scalar{0};
    for (const auto& tet : tets6) {
      std::array<Tetrahedron, 4> subTets;
      decompose(tet, subTets);

      overlapCalcTets24 +=
          overlap_volume(s, std::begin(subTets), std::end(subTets));
    }

    const auto overlapCalcWedges =
        overlap_volume(s, std::begin(wedges), std::end(wedges));

    CHECK(overlapCalcHex == Approx(overlapCalcWedges).epsilon(epsilon));
    CHECK(overlapCalcHex == Approx(overlapCalcTets5).epsilon(epsilon));
    CHECK(overlapCalcHex == Approx(overlapCalcTets6).epsilon(epsilon));
    CHECK(overlapCalcHex == Approx(overlapCalcTets24).epsilon(epsilon));
    CHECK(overlapCalcTets5 == Approx(overlapCalcTets6).epsilon(epsilon));
    CHECK(overlapCalcTets5 == Approx(overlapCalcTets24).epsilon(epsilon));
    CHECK(overlapCalcTets6 == Approx(overlapCalcTets24).epsilon(epsilon));
  });
}

inline void validate_overlap_area(const Sphere& s, const Hexahedron& hex,
                                  const Scalar epsilon) {
  CHECK_NOTHROW({
    std::array<Tetrahedron, 5> tets5;
    std::array<Tetrahedron, 6> tets6;
    std::array<Wedge, 2> wedges;

    decompose(hex, tets5);
    decompose(hex, tets6);
    decompose(hex, wedges);

    const auto areaCalcHex = overlap_area(s, hex)[0];

    auto areaCalcWedges = Scalar{0};
    auto areaCalcTets5 = Scalar{0};
    auto areaCalcTets6 = Scalar{0};
    auto areaCalcTets24 = Scalar{0};

    for (const auto& tet : tets5) {
      areaCalcTets5 += overlap_area(s, tet)[0];
    }

    for (const auto& tet : tets6) {
      std::array<Tetrahedron, 4> subTets;
      decompose(tet, subTets);

      for (const auto& subTet : subTets) {
        areaCalcTets24 += overlap_area(s, subTet)[0];
      }

      areaCalcTets6 += overlap_area(s, tet)[0];
    }

    for (const auto& wedge : wedges) {
      areaCalcWedges += overlap_area(s, wedge)[0];
    }

    CHECK(areaCalcHex == Approx(areaCalcWedges).epsilon(epsilon));
    CHECK(areaCalcHex == Approx(areaCalcTets5).epsilon(epsilon));
    CHECK(areaCalcHex == Approx(areaCalcTets6).epsilon(epsilon));
    CHECK(areaCalcHex == Approx(areaCalcTets24).epsilon(epsilon));
    CHECK(areaCalcTets5 == Approx(areaCalcTets6).epsilon(epsilon));
    CHECK(areaCalcTets5 == Approx(areaCalcTets24).epsilon(epsilon));
    CHECK(areaCalcTets6 == Approx(areaCalcTets24).epsilon(epsilon));
  });
}

}  // namespace overlap

#endif  // OVERLAP_TEST_COMMON_HPP
