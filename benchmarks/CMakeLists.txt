if(NOT TARGET doctest::doctest)
  include(FetchContent)

  fetchcontent_declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.11
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
    GIT_TAG v4.3.11
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL
  )

  fetchcontent_makeavailable(nanobench)
endif()

# list of benchmarks
set(_benchmarks details hex_overlap tet_overlap)

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
