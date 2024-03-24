if(NOT SKBUILD)
  # CMake 3.18 should support this, but there is a bug w.r.t. multiple queries
  # with different components
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.22.0")
    set(_python_components Interpreter Development.Module)
  else()
    set(_python_components Interpreter Development)
  endif()

  find_package(
    Python 3.5
    COMPONENTS ${_python_components}
    REQUIRED
  )
endif()

# add pybind11 shipped as a third-party dependency
add_subdirectory(
  "${PROJECT_SOURCE_DIR}/third_party/pybind11" "third_party/pybind11"
  EXCLUDE_FROM_ALL
)

set(package_name "overlap")

pybind11_add_module(${package_name}_python MODULE ${package_name}.cpp)

set_target_properties(
  ${package_name}_python PROPERTIES OUTPUT_NAME ${package_name}
)

target_link_libraries(
  ${package_name}_python PRIVATE ${package_name}::headers Eigen3::Eigen
)

if(NOT SKBUILD)
  include(GNUInstallDirs)

  string(CONCAT python_install_dir "${CMAKE_INSTALL_LIBDIR}"
                "/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}"
                "/site-packages/${package_name}"
  )

  install(TARGETS ${package_name}_python
          LIBRARY DESTINATION ${python_install_dir} COMPONENT python
  )
else()
  install(TARGETS ${package_name}_python LIBRARY DESTINATION .)
endif()
