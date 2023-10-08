/*!
 * Exact calculation of the overlap volume of spheres and mesh elements.
 * http://dx.doi.org/10.1016/j.jcp.2016.02.003
 *
 * Copyright (C) 2023 Severin Strobl <severin.strobl@fau.de>
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
      .minEpochIterations(25000)
      .run(name, std::forward<F>(func))
      .render(ankerl::nanobench::templates::pyperf(), log);
}
