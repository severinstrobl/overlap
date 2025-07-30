// Copyright (C) 2023 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include <chrono>
#include <fstream>
#include <string>
#include <utility>

#include <nanobench.h>

template<typename F>
inline auto create_benchmark(const std::string& name, F&& func)
    -> ankerl::nanobench::Bench {
  auto log = std::ofstream{name + ".json"};
  return ankerl::nanobench::Bench()
      .title(name)
      .minEpochIterations(25'000)
      .run(name, std::forward<F>(func))
      .render(ankerl::nanobench::templates::pyperf(), log);
}
