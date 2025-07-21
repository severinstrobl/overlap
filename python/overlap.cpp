// Copyright (C) 2021-2025 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#include <array>
#include <cctype>
#include <locale>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeindex>

#include <pybind11/cast.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "overlap/overlap.hpp"

namespace py = pybind11;

using namespace overlap;

namespace {

template<typename... Args>
using overload_cast_ = py::detail::overload_cast_impl<Args...>;

using namespace py::literals;

template<typename Element,
         typename = std::enable_if_t<detail::is_element_v<Element>>>
void create_bindings(py::module& m) {
  static const auto element_names = std::map<std::type_index, std::string>{
      {std::type_index(typeid(Tetrahedron)), "Tetrahedron"},
      {std::type_index(typeid(Wedge)), "Wedge"},
      {std::type_index(typeid(Hexahedron)), "Hexahedron"}};

  const auto& name = element_names.at(std::type_index(typeid(Element)));
  const auto name_lower =
      std::tolower(name.front(), std::locale{}) + name.substr(1);

  static constexpr auto num_vertices = detail::num_vertices<Element>();

  py::class_<Element>(m, name.c_str())
      .def(py::init<std::array<Vector, num_vertices>>(), py::arg("vertices"))
      .def(py::init([](const py::array_t<double>& array) {
             const auto proxy = array.unchecked<2>();
             if (proxy.shape(0) != num_vertices || proxy.shape(1) != 3) {
               throw std::invalid_argument{
                   "invalid shape for vertex list, must be (" +
                   std::to_string(num_vertices) + ", 3)"};
             }

             auto vertices = std::array<Vector, num_vertices>{};
             for (auto v = 0U; v < proxy.shape(0); ++v) {
               vertices[v] = Vector{proxy(v, 0), proxy(v, 1), proxy(v, 2)};
             }

             return Element{std::move(vertices)};
           }),
           py::arg("vertices"))
      .def_readonly("vertices", &Element::vertices,
                    "Return the vertices of the element.")
      .def_readonly("center", &Element::center,
                    "Return the center point of the element.")
      .def_readonly("volume", &Element::volume,
                    "Return the volume of the element.")
      .def_property_readonly(
          "surface_area",
          [](const Element& elem) { return elem.surface_area(); },
          "Return the surface area of the element.");

  m.def(
      "overlap_volume",
      overload_cast_<const Sphere&, const Element&>()(&overlap_volume<Element>),
      "sphere"_a, py::arg{name_lower.c_str()},
      ("Calculate the overlap volume of a sphere and a " + name_lower + ".")
          .c_str());

  m.def("overlap_area",
        overload_cast_<const Sphere&, const Element&>()(&overlap_area<Element>),
        "sphere"_a, py::arg{name_lower.c_str()},
        ("Calculate the overlap area of a sphere and a " + name_lower + ".")
            .c_str());
}

}  // namespace

PYBIND11_MODULE(_overlap, m) {
  py::class_<Sphere>(m, "Sphere")
      .def(py::init<Vector, Scalar>(), py::arg("center"), py::arg("radius"))
      .def_readonly("center", &Sphere::center,
                    "Return the center point of the sphere.")
      .def_readonly("radius", &Sphere::radius,
                    "Return the radius of the sphere.")
      .def_readonly("volume", &Sphere::volume,
                    "Return the volume of the sphere.")
      .def_property_readonly(
          "surface_area", [](const Sphere& s) { return s.surface_area(); },
          "Return the surface area of the sphere.");

  create_bindings<Tetrahedron>(m);
  create_bindings<Wedge>(m);
  create_bindings<Hexahedron>(m);
}
