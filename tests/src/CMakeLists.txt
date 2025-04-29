# handle configuration for generation of coverage information as library
add_library(coverage INTERFACE)

if(OVERLAP_WITH_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "(GNU|Clang)")
    if(NOT CMAKE_BUILD_TYPE MATCHES "(Debug|RelWithDebInfo)")
      message(
        FATAL_ERROR
          "Generation of coverage information requires debug information, set "
          "CMAKE_BUILD_TYPE=Debug or CMAKE_BUILD_TYPE=RelWithDebInfo."
      )
    endif()

    message(STATUS "Enabling generation of coverage reports for C++ unit tests")

    set(_coverage_compile_flags "--coverage")

    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-fprofile-abs-path" have_fprofile_abs_path)
    if(have_fprofile_abs_path)
      list(APPEND _coverage_compile_flags "-fprofile-abs-path")
    endif()

    target_compile_options(coverage INTERFACE ${_coverage_compile_flags})
    target_link_options(coverage INTERFACE "--coverage")
  else()
    message(FATAL_ERROR "Generation of coverage information not supported for "
                        ${CMAKE_CXX_COMPILER_ID} " compiler"
    )
  endif()
endif()

if(NOT TARGET doctest::doctest)
  include(FetchContent)

  set(_doctest_patch "${CMAKE_CURRENT_BINARY_DIR}/doctest_include.patch")
  file(
    DOWNLOAD
    https://github.com/severinstrobl/doctest/compare/master...patchset.patch
    ${_doctest_patch}
    EXPECTED_HASH
      SHA256=31b3e5004a6420db0799d39c3fa8199d3cd23b388d94d21c324dfbd01ccce93b
  )

  fetchcontent_declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG 1da23a3e8119ec5cce4f9388e91b065e20bf06f5 # v2.4.12
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL
    PATCH_COMMAND patch -p1 --forward < ${_doctest_patch} || true
  )

  fetchcontent_makeavailable(doctest)
  list(APPEND CMAKE_MODULE_PATH ${doctest_SOURCE_DIR}/scripts/cmake)
endif()

include(CTest)
include(doctest)

# helper function to add an individual test
function(overlap_add_test test_name test_source)
  add_executable(${test_name} ${test_source})
  target_compile_features(${test_name} PRIVATE cxx_std_17)
  target_link_libraries(${test_name} PRIVATE doctest::doctest coverage)

  set_target_properties(${test_name} PROPERTIES FOLDER tests)

  target_compile_definitions(
    ${test_name} PRIVATE DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
  )

  doctest_discover_tests(
    ${test_name}
    REPORTER
    junit
    OUTPUT_DIR
    ${CMAKE_BINARY_DIR}/
    OUTPUT_PREFIX
    doctest_testing_
    OUTPUT_SUFFIX
    .xml
    ADD_LABELS
    1
  )
endfunction()

# list of unit tests
set(_unit_tests
    clamp
    contains
    decompose_elements
    detect_non_planar_faces
    double_precision
    elements
    general_wedge
    line_sphere_intersection
    normal_newell
    normalize_element
    polygon
    regularized_wedge
    regularized_wedge_area
    sphere
    sphere_element_area_edgecases
    sphere_element_overlap
    sphere_element_overlap_edgecases
    sphere_hex_area
    sphere_tet_area
    sphere_tet_overlap_edgecases
    unit_sphere_intersections
)

list(TRANSFORM _unit_tests PREPEND "test_")

# check for feenableexcept
include(CheckIncludeFile)
include(CheckSymbolExists)
include(CMakePushCheckState)

cmake_push_check_state(RESET)
check_include_file("fenv.h" _have_fenv_h)
if(_have_fenv_h)
  set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
  if(UNIX)
    set(CMAKE_REQUIRED_LIBRARIES m)
  endif()
  check_symbol_exists(feenableexcept "fenv.h" OVERLAP_HAVE_FEENABLEEXCEPT)
  set(CMAKE_REQUIRED_LIBRARIES)
endif()
cmake_reset_check_state()

# activate additional compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "(GNU|Clang)" AND CMAKE_BUILD_TYPE MATCHES
                                                   "(Debug|RelWithDebInfo)"
)
  string(APPEND CMAKE_CXX_FLAGS " -Wall -Werror -Wextra")
endif()

# register the individual unit tests
foreach(unit_test ${_unit_tests})
  overlap_add_test(${unit_test} "${unit_test}.cpp")
  target_link_libraries(${unit_test} PRIVATE overlap::headers Eigen3::Eigen)

  if(OVERLAP_HAVE_FEENABLEEXCEPT)
    target_compile_definitions(
      ${unit_test} PRIVATE -DOVERLAP_HAVE_FEENABLEEXCEPT
    )
  endif()
endforeach()
