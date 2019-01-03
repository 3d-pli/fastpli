CMAKE_MINIMUM_REQUIRED(VERSION 3.2)

PROJECT( fastPLI )
SET(CMAKE_COLOR_MAKEFILE ON)
# SET(CMAKE_VERBOSE_MAKEFILE ON)

##########################################################################
# Buildconfiguration
##########################################################################
IF( NOT CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE Debug )
ENDIF()

SET(flags "-std=c++17 -O3 -march=native")
# SET(flags "-std=c++17 -g")

IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(warnings "-Wall -Wpedantic -Wextra -Wfatal-errors")
ENDIF()

SET(CMAKE_CXX_FLAGS "${warnings} ${flags}")

##########################################################################
# Libraries
##########################################################################
ADD_SUBDIRECTORY(lib/pybind11)

# OpenMP
FIND_PACKAGE( OpenMP )
IF ( OPENMP_FOUND )
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ELSE ( OPENMP_FOUND )
    MESSAGE("OPENMP NOT FOUND")
ENDIF ( OPENMP_FOUND )

# OpenGL
FIND_PACKAGE( OpenGL )
FIND_PACKAGE( GLUT )

##########################################################################
# Sources
##########################################################################
INCLUDE_DIRECTORIES(src)
ADD_SUBDIRECTORY("${PROJECT_SOURCE_DIR}/src/model")
ADD_SUBDIRECTORY("${PROJECT_SOURCE_DIR}/src/objects")
