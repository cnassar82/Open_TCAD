# SUPREM geometry/mesh core (M0, suprem4/doc/plan_m0_geometry_core.md)
cmake_minimum_required(VERSION 3.16)
project(suprem_core C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# Mmg v5.8.0 (LGPL-3), 2D library only, static, no optional dependencies.
set(BUILD MMG2D CACHE STRING "" FORCE)
set(LIBMMG2D_STATIC ON CACHE BOOL "" FORCE)
set(LIBMMG2D_SHARED OFF CACHE BOOL "" FORCE)
set(USE_SCOTCH OFF CACHE BOOL "" FORCE)
set(USE_ELAS OFF CACHE BOOL "" FORCE)
set(USE_VTK OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
add_subdirectory(third_party/mmg EXCLUDE_FROM_ALL)

set(CLIPPER2_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src/third_party/clipper2)
add_library(suprem_core STATIC
  src/str_file.cpp
  src/mesh.cpp
  src/brep.cpp
  src/locate.cpp
  src/remesh.cpp
  src/transfer.cpp
  src/evolve.cpp
  src/capi.cpp
  ${CLIPPER2_DIR}/src/clipper.engine.cpp
  ${CLIPPER2_DIR}/src/clipper.offset.cpp
  ${CLIPPER2_DIR}/src/clipper.rectclip.cpp)
target_include_directories(suprem_core PUBLIC include
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cdt/include
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/third_party/clipper2/include)
target_include_directories(suprem_core PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/third_party/mmg/include)
target_link_libraries(suprem_core PUBLIC libmmg2d_a)
target_compile_options(suprem_core PRIVATE -Wall -Wextra)

add_executable(core_str tools/core_str.cpp)
target_link_libraries(core_str suprem_core)

add_executable(core_remesh tools/core_remesh.cpp)
target_link_libraries(core_remesh suprem_core)

add_executable(core_evolve_test tools/core_evolve_test.cpp)
target_link_libraries(core_evolve_test suprem_core)
