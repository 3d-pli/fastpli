CMAKE_MINIMUM_REQUIRED(VERSION 3.2)

PROJECT( fastpli LANGUAGES CXX)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

##########################################################################
# Buildconfiguration
##########################################################################
if( NOT CMAKE_BUILD_TYPE )
   set(CMAKE_BUILD_TYPE "Debug")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   set(warnings "-Wall -Wpedantic -Wextra -Wfatal-errors")
   set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${warnings})
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
endif()

##########################################################################
# BUILD VERSION FILES
##########################################################################
configure_file(
   ${CMAKE_SOURCE_DIR}/cmake/version.cmake
   ${CMAKE_CURRENT_BINARY_DIR}/configure_version.cmake
)
add_custom_target(
   build_version
   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/configure_version.cmake
)

##########################################################################
# Libraries
##########################################################################
# Pybind11
FIND_PACKAGE(PythonInterp 3 REQUIRED)
FIND_PACKAGE(PythonLibs   3 REQUIRED)
ADD_SUBDIRECTORY(lib/pybind11)

# OpenMP
FIND_PACKAGE( OpenMP )
IF ( OPENMP_FOUND )
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ELSE ( OPENMP_FOUND )
    MESSAGE("OPENMP NOT FOUND")
ENDIF ( OPENMP_FOUND )

# MPI
FIND_PACKAGE( MPI REQUIRED )
INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_DIRS})

# OpenGL
set(OpenGL_GL_PREFERENCE "GLVND")
FIND_PACKAGE( OpenGL )
FIND_PACKAGE( GLUT )

##########################################################################
# Sources
##########################################################################
INCLUDE_DIRECTORIES( AFTER src)
ADD_SUBDIRECTORY("${PROJECT_SOURCE_DIR}/src/model/solver")
ADD_SUBDIRECTORY("${PROJECT_SOURCE_DIR}/src/objects")
ADD_SUBDIRECTORY("${PROJECT_SOURCE_DIR}/src/simulation")
