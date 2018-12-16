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
    SET(warnings "-Wall -Wpedantic -Wextra -Werror -Wfatal-errors")
ENDIF()

SET(CMAKE_CXX_FLAGS "${warnings} ${flags}")

##########################################################################
# Libraries
##########################################################################
ADD_SUBDIRECTORY(lib/pybind11)

##########################################################################
# Sources
##########################################################################
INCLUDE_DIRECTORIES(src/include)
ADD_SUBDIRECTORY("${PROJECT_SOURCE_DIR}/src/model/solver")
ADD_SUBDIRECTORY("${PROJECT_SOURCE_DIR}/src/objects")
