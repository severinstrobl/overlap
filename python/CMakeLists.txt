if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.18.0")
  set(python_components Interpreter Development.Module)
else()
  set(python_components Interpreter Development)
endif()

# get FindPython working with scikit-build
if(SKBUILD)
  set(Python_EXECUTABLE "${PYTHON_EXECUTABLE}")
  set(Python_INCLUDE_DIR "${PYTHON_INCLUDE_DIR}")
  unset(PYTHON_EXECUTABLE)
  unset(PYTHON_INCLUDE_DIR)
  unset(PYTHON_LIBRARY)
  unset(PYTHON_VERSION_STRING)
endif()

find_package(
  Python 3.6
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

if(CMAKE_PROJECT_VERSION)
target_compile_definitions(_overlap
                           PRIVATE VERSION_INFO="${CMAKE_PROJECT_VERSION}")
endif()

install(TARGETS _overlap DESTINATION .)
