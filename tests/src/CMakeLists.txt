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

# FIXME: switch to Catch2

# build GoogleTest shipped as a third-party dependency

# avoid linker issues on Windows when mixing shared and static libraries
set(gtest_force_shared_crt
    ON
    CACHE BOOL "" FORCE
)

add_subdirectory(
  "${PROJECT_SOURCE_DIR}/third_party/googletest" "third_party/googletest"
  EXCLUDE_FROM_ALL
)

mark_as_advanced(
  BUILD_GMOCK
  BUILD_GTEST
  BUILD_SHARED_LIBS
  gmock_build_tests
  gtest_build_samples
  gtest_build_tests
  gtest_disable_pthreads
  gtest_force_shared_crt
  gtest_hide_internal_symbols
)

set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)
set_target_properties(gmock PROPERTIES FOLDER extern)
set_target_properties(gmock_main PROPERTIES FOLDER extern)

# load Google Test module and set up helper function for creating tests
include(GoogleTest)

# helper function to dd a test using the Google Test framework
function(overlap_add_test test_name test_source)
  add_executable(${test_name} ${test_source})
  target_link_libraries(${test_name} PRIVATE gtest gtest_main coverage)

  gtest_discover_tests(
    ${test_name}
    WORKING_DIRECTORY ${PROJECT_DIR}
    PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
  )

  set_target_properties(${test_name} PROPERTIES FOLDER tests)
endfunction()

# list of unit tests
set(unit_tests
    clamp
    decompose_elements
    elements
    normal_newell
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
)

list(TRANSFORM unit_tests PREPEND "test_")

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

# register the individual unit tests
foreach(unit_test ${unit_tests})
  overlap_add_test(${unit_test} "${unit_test}.cpp")
  target_link_libraries(${unit_test} PRIVATE overlap::headers Eigen3::Eigen)

  if(OVERLAP_HAVE_FEENABLEEXCEPT)
    target_compile_definitions(
      ${unit_test} PRIVATE -DOVERLAP_HAVE_FEENABLEEXCEPT
    )
  endif()
endforeach()
