if(NOT TARGET doctest::doctest)
  include(FetchContent)

  fetchcontent_declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG 415855d33e88146425c2d8db1825797e7b1a4c85 # v2.5.1
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL
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
set(_benchmarks details element_construction hex_overlap_volume
                tet_overlap_volume
)

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
