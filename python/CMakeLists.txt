# CMake 3.18 should support this, but there is a bug w.r.t. multiple queries
# with different components
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.24.0")
  set(python_components Interpreter Development.Module)
else()
  set(python_components Interpreter Development)
endif()

find_package(
  Python 3.5
  COMPONENTS ${python_components}
  REQUIRED
)

# build pybind11 shipped as a third-party dependency
add_subdirectory(
  "${PROJECT_SOURCE_DIR}/third_party/pybind11" "third_party/pybind11"
  EXCLUDE_FROM_ALL
)

pybind11_add_module(_overlap MODULE overlap.cpp)

target_link_libraries(_overlap PRIVATE Eigen3::Eigen)

install(TARGETS _overlap DESTINATION .)
