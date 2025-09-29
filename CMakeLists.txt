cmake_minimum_required(VERSION 3.10)

project(overset_mesh_solver CXX)

# Require at least C++11 for GCC 6.x compatibility
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Sources (new folder structure)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SOURCES
    ${SRC_DIR}/main.cpp
    ${SRC_DIR}/mesh/adt.cpp
    ${SRC_DIR}/mesh/datastructure.cpp
    ${SRC_DIR}/solver/solver.cpp
    ${SRC_DIR}/util/utilities.cpp
    ${SRC_DIR}/io/output.cpp
)

add_executable(overset ${SOURCES})

target_include_directories(overset PRIVATE 
  ${SRC_DIR}
  ${SRC_DIR}/mesh
  ${SRC_DIR}/solver
  ${SRC_DIR}/util
  ${SRC_DIR}/io
)

# Helpful warnings for GCC/Clang
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(overset PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Strong optimization flags in Release
if (CMAKE_BUILD_TYPE STREQUAL "Release")
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(overset PRIVATE -O3 -DNDEBUG)
  endif()
endif()
