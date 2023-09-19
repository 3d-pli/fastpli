cmake_minimum_required(VERSION 3.9)

project(fastpli LANGUAGES CXX)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# ##############################################################################
# Buildconfiguration
# ##############################################################################
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# warning IF not release build
if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  set(warn "BUILD TYPE: ${CMAKE_BUILD_TYPE}")
  if(NOT WIN32)
    string(ASCII 27 Esc)
    set(warn "${Esc}[1;33m${warn}${Esc}[m")
  endif()
  message(WARNING "${warn}")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# check architecture
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH)
if(COMPILER_SUPPORTS_MARCH)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
else()
  unset(COMPILER_SUPPORTS_MARCH)
endif()

if(NOT COMPILER_SUPPORTS_MARCH)
  check_cxx_compiler_flag("-mcpu=native" COMPILER_SUPPORTS_MCPU)
  if(COMPILER_SUPPORTS_MCPU)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native")
  endif()
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}"
                                                STREQUAL "Clang")
  set(warnings "-Wall -Wpedantic -Wextra -Wfatal-errors -Wdeprecated -Werror")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warnings}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -O3 -g3")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb -g3")
else()
  message(WARNING "Compiler not supported yet, deactivating optimization flags")
endif()

list(APPEND FASTPLI_INFO_COMPILER "${CMAKE_CXX_COMPILER_ID}
                                    ${CMAKE_CXX_COMPILER}
                                    ${CMAKE_CXX_COMPILER_VERSION}")

# ##############################################################################
# Git Submodules
# ##############################################################################
find_package(Git REQUIRED)
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/.git")
  message(
    FATAL_ERROR
      "fastPLI compromised. .git not found. Please clone the Project again.")
else()
  message(STATUS "Submodule update")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE GIT_SUBMOD_RESULT)
  if(NOT GIT_SUBMOD_RESULT EQUAL "0")
    message(
      FATAL_ERROR string ("git submodule update --init failed with"
                          "${GIT_SUBMOD_RESULT}, please checkout submodules"))
  endif()
endif()

# ##############################################################################
# Libraries
# ##############################################################################
# Python

# SET(PYTHON3_VERSION "" CACHE STRING "")

message("Python3_VERSION: ${Python3_VERSION}")

if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
  find_package(PythonInterp EXACT ${Python3_VERSION} REQUIRED)
  set(Python3_VERSION
      "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}.${PYTHON_VERSION_PATCH}")
  find_package(PythonLibs EXACT ${Python3_VERSION} REQUIRED)
else()
  find_package(
    Python3 EXACT ${Python3_VERSION}
    COMPONENTS Interpreter Development
    REQUIRED)
endif()
list(APPEND FASTPLI_INFO_LIBRARIES "Python-v${Python3_VERSION}")

# Pybind11
add_subdirectory(lib/pybind11)
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsized-deallocation")
endif()
include_directories(${Python3_INCLUDE_DIRS})
list(APPEND FASTPLI_INFO_LIBRARIES "Pybind-v${PYBIND11_PYTHON_VERSION}")

# OpenMP
find_package(OpenMP)
if(OPENMP_CXX_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else(OPENMP_CXX_FOUND)
  message("OPENMP NOT FOUND")
endif(OPENMP_CXX_FOUND)
list(APPEND FASTPLI_INFO_LIBRARIES "OpenMP_CXX-v${OpenMP_CXX_VERSION}")

# MPI
find_package(MPI REQUIRED)
list(APPEND FASTPLI_INFO_LIBRARIES "MPI_CXX-v${MPI_CXX_VERSION}")

# OpenGL
set(OpenGL_GL_PREFERENCE "GLVND")
find_package(OpenGL)
find_package(GLUT)
list(APPEND FASTPLI_INFO_LIBRARIES "OpenGL-v${OpenGL_VERSION}")
list(APPEND FASTPLI_INFO_LIBRARIES "GLUT-v${GLUT_VERSION}")

# ##############################################################################
# BUILD VERSION FILES
# ##############################################################################
configure_file(${CMAKE_SOURCE_DIR}/cmake/init.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/configure_init.cmake)

add_custom_target(
  init_files ALL COMMAND ${CMAKE_COMMAND} -P
                         ${CMAKE_CURRENT_BINARY_DIR}/configure_init.cmake)

# ##############################################################################
# Sources
# ##############################################################################
include_directories(AFTER src)

set(PROJECT_LIB_DIR "${PROJECT_SOURCE_DIR}/src/fastpli")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/model/solver")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/objects")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/simulation")
