# build GoogleTest shipped as a third-party dependency.

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
  set_target_properties(${test_name} PROPERTIES FOLDER test)
endfunction()

# list of unit tests
set(unit_tests
    decompose_elements
    detail_clamp
    detail_regularizedWedge
    detail_regularizedWedgeArea
    elements
    normal_newell
    polygon
    sphere
    sphere_element_area_edgecases
    sphere_element_overlap
    sphere_element_overlap_edgecases
    sphere_hex_area
    sphere_tet_area
    sphere_tet_overlap_edgecases
)

# check for feenableexcept
include(CheckSymbolExists)
check_symbol_exists(feenableexcept "fenv.h" HAVE_FEENABLEEXCEPT)

# register the individual unit tests
foreach(unit_test ${unit_tests})
  overlap_add_test(${unit_test} "${unit_test}.cpp")
  target_link_libraries(${unit_test} PUBLIC Eigen3::Eigen)

  if(HAVE_FEENABLEEXCEPT)
    target_compile_definitions(${unit_test} PRIVATE -DHAVE_FEENABLEEXCEPT)
  endif()
endforeach()
