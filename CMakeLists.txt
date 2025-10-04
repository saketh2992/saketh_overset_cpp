cmake_minimum_required(VERSION 3.10)

project(overset_mesh_solver CXX)

# Require at least C++11 for GCC 6.x compatibility
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Sources (new folder structure)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Common source files (excluding main)
set(COMMON_SOURCES
    ${SRC_DIR}/mesh/adt.cpp
    ${SRC_DIR}/mesh/datastructure.cpp
    ${SRC_DIR}/mesh/datastructure_rectangular.cpp
    ${SRC_DIR}/mesh/datastructure_circular.cpp
    ${SRC_DIR}/solver/solver.cpp
    ${SRC_DIR}/solver/solver_circular.cpp
    ${SRC_DIR}/util/utilities.cpp
    ${SRC_DIR}/io/output.cpp
)

# Original executable - DISABLED (solver not yet updated for inheritance)
# set(SOURCES ${SRC_DIR}/main.cpp ${COMMON_SOURCES})
# add_executable(overset ${SOURCES})

# New executable with inheritance-based implementation
set(SOURCES_NEW ${SRC_DIR}/main_new.cpp ${COMMON_SOURCES})
add_executable(overset_new ${SOURCES_NEW})

target_include_directories(overset_new PRIVATE 
  ${SRC_DIR}
  ${SRC_DIR}/mesh
  ${SRC_DIR}/solver
  ${SRC_DIR}/util
  ${SRC_DIR}/io
)

# Cylinder flow executable (rectangular background + circular component)
set(SOURCES_CYLINDER ${SRC_DIR}/main_cylinder_flow.cpp ${COMMON_SOURCES})
add_executable(overset_cylinder ${SOURCES_CYLINDER})

target_include_directories(overset_cylinder PRIVATE 
  ${SRC_DIR}
  ${SRC_DIR}/mesh
  ${SRC_DIR}/solver
  ${SRC_DIR}/util
  ${SRC_DIR}/io
)

# Helpful warnings for GCC/Clang
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # target_compile_options(overset PRIVATE -Wall -Wextra -Wpedantic)  # DISABLED
  target_compile_options(overset_new PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(overset_cylinder PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Strong optimization flags in Release
if (CMAKE_BUILD_TYPE STREQUAL "Release")
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # target_compile_options(overset PRIVATE -O3 -DNDEBUG)  # DISABLED
    target_compile_options(overset_new PRIVATE -O3 -DNDEBUG)
    target_compile_options(overset_cylinder PRIVATE -O3 -DNDEBUG)
  endif()
endif()

# ============================================================================
# OPTIONAL: TIOGA Overset Grid Assembly Library
# ============================================================================
# Find MPI (required for TIOGA)
find_package(MPI QUIET)

# Set TIOGA installation directory
get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
set(TIOGA_DIR "${PARENT_DIR}/tioga_install" CACHE PATH "TIOGA installation directory")

# Try to find TIOGA library and headers
find_library(TIOGA_LIBRARY 
    NAMES tioga 
    PATHS ${TIOGA_DIR}/lib 
          $ENV{HOME}/tioga_install/lib
    NO_DEFAULT_PATH
)

find_path(TIOGA_INCLUDE_DIR 
    NAMES tioga.h 
    PATHS ${TIOGA_DIR}/include
          $ENV{HOME}/tioga_install/include
    NO_DEFAULT_PATH
)

# Optionally link TIOGA if found
if(TIOGA_LIBRARY AND TIOGA_INCLUDE_DIR AND MPI_FOUND)
    message(STATUS "═══════════════════════════════════════════════════")
    message(STATUS "TIOGA Library Found!")
    message(STATUS "  Library: ${TIOGA_LIBRARY}")
    message(STATUS "  Include: ${TIOGA_INCLUDE_DIR}")
    message(STATUS "  MPI:     ${MPI_CXX_LIBRARIES}")
    message(STATUS "═══════════════════════════════════════════════════")
    
    # Define preprocessor flag for conditional compilation
    add_compile_definitions(HAS_TIOGA)
    
    # Add TIOGA and MPI include directories
    # target_include_directories(overset PRIVATE  ${TIOGA_INCLUDE_DIR} ${MPI_CXX_INCLUDE_DIRS})  # DISABLED
    target_include_directories(overset_new PRIVATE 
        ${TIOGA_INCLUDE_DIR}
        ${MPI_CXX_INCLUDE_DIRS}
    )
    target_include_directories(overset_cylinder PRIVATE 
        ${TIOGA_INCLUDE_DIR}
        ${MPI_CXX_INCLUDE_DIRS}
    )
    
    # Link TIOGA and MPI libraries
    # target_link_libraries(overset ${TIOGA_LIBRARY} ${MPI_CXX_LIBRARIES})  # DISABLED
    target_link_libraries(overset_new 
        ${TIOGA_LIBRARY}
        ${MPI_CXX_LIBRARIES}
    )
    target_link_libraries(overset_cylinder 
        ${TIOGA_LIBRARY}
        ${MPI_CXX_LIBRARIES}
    )
    
    # Add MPI compile flags if needed
    if(MPI_CXX_COMPILE_FLAGS)
        # target_compile_options(overset PRIVATE ${MPI_CXX_COMPILE_FLAGS})  # DISABLED
        target_compile_options(overset_new PRIVATE ${MPI_CXX_COMPILE_FLAGS})
        target_compile_options(overset_cylinder PRIVATE ${MPI_CXX_COMPILE_FLAGS})
    endif()
    
    # Add MPI link flags if needed
    if(MPI_CXX_LINK_FLAGS)
        # target_link_options(overset PRIVATE ${MPI_CXX_LINK_FLAGS})  # DISABLED
        target_link_options(overset_new PRIVATE ${MPI_CXX_LINK_FLAGS})
        target_link_options(overset_cylinder PRIVATE ${MPI_CXX_LINK_FLAGS})
    endif()
    
else()
    if(NOT MPI_FOUND)
        message(STATUS "MPI not found - building without TIOGA support")
    elseif(NOT TIOGA_LIBRARY)
        message(STATUS "TIOGA library not found at: ${TIOGA_DIR}/lib")
        message(STATUS "Building without overset grid assembly support")
    elseif(NOT TIOGA_INCLUDE_DIR)
        message(STATUS "TIOGA headers not found at: ${TIOGA_DIR}/include")
        message(STATUS "Building without overset grid assembly support")
    endif()
    message(STATUS "To install TIOGA, run: ./install_tioga.sh")
endif()
