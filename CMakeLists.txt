cmake_minimum_required(VERSION 3.10)

project(overset_mesh_solver CXX)

# Require at least C++11 for GCC 6.x compatibility
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Sources
set(SOURCES
    main.cpp
    adt.cpp
    datastructure.cpp
    solver.cpp
    utilities.cpp
    output.cpp
)

add_executable(overset ${SOURCES})

target_include_directories(overset PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

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
