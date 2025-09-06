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

if(NOT TARGET nanobench)
  include(FetchContent)

  fetchcontent_declare(
    nanobench
    GIT_REPOSITORY https://github.com/martinus/nanobench.git
    GIT_TAG a5a50c2b33eea2ff1fcb355cacdface43eb42b25 # v4.3.11
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL
  )

  fetchcontent_makeavailable(nanobench)
endif()

# list of benchmarks
set(_benchmarks details hex_overlap_volume tet_overlap_volume)

# register the individual benchmarks
foreach(benchmark ${_benchmarks})
  add_executable(${benchmark} ${benchmark}.cpp)
  target_link_libraries(
    ${benchmark} PRIVATE overlap::headers Eigen3::Eigen nanobench
                         doctest::doctest
  )

  target_compile_definitions(
    ${benchmark} PRIVATE DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
  )

  add_test(NAME ${benchmark} COMMAND $<TARGET_FILE:${benchmark}>)

  set_tests_properties(${benchmark} PROPERTIES LABELS benchmark)
endforeach()
