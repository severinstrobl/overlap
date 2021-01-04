# Build pybind11 shipped as a third-party dependency.
add_subdirectory("${PROJECT_SOURCE_DIR}/third_party/pybind11"
                 "third_party/pybind11"
                 EXCLUDE_FROM_ALL)

pybind11_add_module(_overlap MODULE overlap.cpp)

target_link_libraries(_overlap PRIVATE Eigen3::Eigen)

# Forward version info from setup.py.
#target_compile_definitions(overlap PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

install(TARGETS _overlap DESTINATION .)