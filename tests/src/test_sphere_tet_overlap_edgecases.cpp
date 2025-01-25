// Copyright (C) 2015-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("SphereTetOverlap") {
  using namespace overlap;

  TEST_CASE("EdgeCases") {
    // clang-format off
		const auto spheres = std::vector<Sphere>{
			{Vector::Zero(), 1.0},
			{{-0.01725, 0, 0}, 1.0},
			{Vector::Zero(), 1.0}
		};

		const auto tetrahedra = std::vector<Tetrahedron>{
			{{{
				{0.0357829, 0, 1.01271},
				{0, 0, 1.01271},
				{0.0356948, 0.0386086, 0.962075},
				{0, 0, 0.962075}
			}}},
			{{{
				{0.9667906976744187, 0, 3.098296812907414e-16},
				{1.002654107311333, 0.0384643285369352, -2.82302142880589e-16},
				{1.002573643410853, 0, 4.131062417209885e-16},
				{1.002573643410853, 0, -0.05063534883720874}
			}}},
			{{{
				{0.28, -0.9599999999999999, -0.02102000000000028},
				{0.2400000000000001, -0.9599999999999999, 0.01898000000000015},
				{0.28, -0.9999999999999999, 0.01898000000000015},
				{0.28, -0.9599999999999999, 0.01898000000000015}
			}}}
		};
    // clang-format on

    REQUIRE(spheres.size() == tetrahedra.size());

    const auto epsilon = std::sqrt(std::numeric_limits<Scalar>::epsilon());
    for (std::size_t idx = 0; idx < spheres.size(); ++idx) {
      std::array<Tetrahedron, 4> tets4;
      decompose(tetrahedra[idx], tets4);

      const auto overlapCalcTet = overlap_volume(spheres[idx], tetrahedra[idx]);
      auto overlapCalcTets4 = Scalar{0};
      for (const auto& tet : tets4) {
        overlapCalcTets4 += overlap_volume(spheres[idx], tet);
      }

      CHECK(overlapCalcTet ==
            Approx(overlapCalcTets4).epsilon(epsilon * spheres[idx].volume));
    }
  }

  TEST_CASE("#96") {
    const auto center = Vector{1.7553357e-06, 4.2232066e-06, 5.8329073e-07};
    const auto sphere = Sphere{center, 20e-9};

    const auto vertices = std::array<Vector, 4>{
        {{1.7503302395906002e-06, 4.2330364312997e-06, 5.961778422123901e-07},
         {1.7438173901207002e-06, 4.222375361573301e-06, 5.9263766042144e-07},
         {1.7394539738699001e-06, 4.2382759184772e-06, 6.009593818316999e-07},
         {1.7544257028301e-06, 4.2350646020068004e-06, 5.840237397166e-07}}};

    const auto tet = Tetrahedron{vertices};

    CHECK(overlap_volume(sphere, tet) >= 0.0);
  }

  TEST_CASE("#104-inside") {
    const auto center = Vector{5.009999999999999e-07, 5.2e-07, 5e-7};
    const auto sphere = Sphere{center, 20e-9};

    const auto vertices = std::array<Vector, 4>{{{5e-7, 1e-6, 5e-7},
                                                 {1e-6, 5e-7, 5e-7},
                                                 {0, 5e-7, 5e-7},
                                                 {5e-7, 5e-7, 0}}};

    const auto tet = Tetrahedron{vertices};

    CHECK(overlap_volume(sphere, tet) == Approx(0.5 * sphere.volume));
  }

  TEST_CASE("#104-outside") {
    const auto center = Vector{5.009999999999999e-07, 5.2e-07, 5e-7};
    const auto sphere = Sphere{center, 20e-9};

    const auto vertices = std::array<Vector, 4>{{{0, 5e-7, 5e-7},
                                                 {5e-7, 0, 5e-7},
                                                 {1e-6, 5e-7, 5e-7},
                                                 {5e-7, 5e-7, 1e-6}}};

    const auto tet = Tetrahedron{vertices};

    CHECK(overlap_volume(sphere, tet) == Approx(0.5 * sphere.volume));
  }
}
