// Copyright (C) 2015-2022 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include "common.hpp"

TEST_SUITE("SphereElementOverlap") {
  using namespace overlap;

  TEST_CASE("EdgeCases") {
    const auto spheres = std::vector<Sphere>{
        {{0, 0, 0}, 1.0},
        {{-1, 3, 9}, 10},
        {{0.5, 0.5, 0.5}, 0.50001 * std::sqrt(2.0)},

        {{-1.534427712524021, -0.6526040637766801, 3.823443102163421},
         5.459817873898927},

        {{-2.291983426015874, -3.495618444307236, 2.067917670011271},
         4.797942866073771},

        {{-0.2174878528692581, -3.076535346840716, 0.53771818665538},
         2.856370661961459},

        {{-0.4599350500370143, 4.782461234807632, -2.6269327661783},
         5.2681424105015},

        {{-0.7611917089641156, -0.8319982272779169, -0.004847761469840783},
         2.103084880441632},

        {{2.992123379449451, -0.4987719594414469, 1.44196971013958},
         4.706537474211725},

        {{15.461110118225834, -8.575216180676412, 14.487981174363448},
         21.97121086612232}};

    const auto epsilon = std::sqrt(std::numeric_limits<Scalar>::epsilon());
    for (const Sphere& sphere : spheres) {
      validate_overlap_volume(sphere, unit_hexahedron(),
                              epsilon * sphere.volume);
    }
  }
}
