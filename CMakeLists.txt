cmake_minimum_required(VERSION 3.20)

project(coldet VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- Library: coldet (static) ----
set(COLDET_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Glob all source + header files under src/
file(GLOB_RECURSE COLDET_SOURCES CONFIGURE_DEPENDS
  "${COLDET_SRC_DIR}/*.cpp"
  "${COLDET_SRC_DIR}/*.cc"
  "${COLDET_SRC_DIR}/*.cxx"
  "${COLDET_SRC_DIR}/*.c"
  "${COLDET_SRC_DIR}/*.hpp"
  "${COLDET_SRC_DIR}/*.hh"
  "${COLDET_SRC_DIR}/*.hxx"
  "${COLDET_SRC_DIR}/*.h"
)

add_library(coldet STATIC ${COLDET_SOURCES})

target_include_directories(coldet
  PUBLIC
    "${COLDET_SRC_DIR}"
)

# ---- Optional dependency: TBB ----
find_package(TBB QUIET)
find_package(Eigen3 REQUIRED)

if(TBB_FOUND)
  target_compile_definitions(coldet PUBLIC USE_TBB=1)
  target_link_libraries(coldet PUBLIC TBB::tbb)
endif()

# ---- Tests (gtest) ----
include(CTest)
if(BUILD_TESTING)
  find_package(GTest REQUIRED)

  add_executable(coldet_unit_tests
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/main.cpp"
  )

  target_link_libraries(coldet_unit_tests
    PRIVATE
      coldet
      GTest::gtest_main
  )

  # If you want USE_TBB visible to the test target even if it doesn't link coldet publicly:
  if(TBB_FOUND)
    target_compile_definitions(coldet_unit_tests PRIVATE USE_TBB=1)
  endif()

endif()

add_subdirectory(obb)
add_subdirectory(benchmark)
add_subdirectory(tests)

