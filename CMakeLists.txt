CMAKE_MINIMUM_REQUIRED(VERSION 3.9)

PROJECT( fastpli LANGUAGES CXX)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

##########################################################################
# Buildconfiguration
##########################################################################
IF( NOT CMAKE_BUILD_TYPE )
   SET(CMAKE_BUILD_TYPE "Release")
ENDIF()

# warning IF not release build
IF( NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release" )
   SET(warn "BUILD TYPE: ${CMAKE_BUILD_TYPE}")
   IF(NOT WIN32)
      string(ASCII 27 Esc)
      SET(warn "${Esc}[1;33m${warn}${Esc}[m")
   ENDIF()
   MESSAGE(WARNING "${warn}")
ENDIF()

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   SET(warnings "-Wall -Wpedantic -Wextra -Wfatal-errors -Wdeprecated -Werror")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warnings}")
   SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
   SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
   SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
ENDIF()

LIST (APPEND FASTPLI_INFO_COMPILER "${CMAKE_CXX_COMPILER_ID}
                                    ${CMAKE_CXX_COMPILER}
                                    ${CMAKE_CXX_COMPILER_VERSION}" )

##########################################################################
# Git Submodules
##########################################################################
FIND_PACKAGE(Git REQUIRED)
IF(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    OPTION(GIT_SUBMODULE "Check submodules during build" ON)
    IF(GIT_SUBMODULE)
        MESSAGE(STATUS "Submodule update")
        EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        IF(NOT GIT_SUBMOD_RESULT EQUAL "0")
            MESSAGE(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        ENDIF()
    ENDIF()
ENDIF()

##########################################################################
# Libraries
##########################################################################
# Python
FIND_PACKAGE(PythonInterp 3.6 REQUIRED)
FIND_PACKAGE(PythonLibs   3.6 REQUIRED)

ADD_SUBDIRECTORY(lib/pybind11)
IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsized-deallocation" )
ENDIF()
LIST(APPEND FASTPLI_INFO_LIBRARIES "PythonInterp-v${PythonInterp_VERSION}")
LIST(APPEND FASTPLI_INFO_LIBRARIES "PythonLibs-v${PythonLibs_VERSION}")

# OpenMP
FIND_PACKAGE( OpenMP )
IF ( OPENMP_CXX_FOUND )
   SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ELSE ( OPENMP_CXX_FOUND )
   MESSAGE("OPENMP NOT FOUND")
ENDIF ( OPENMP_CXX_FOUND )
LIST(APPEND FASTPLI_INFO_LIBRARIES "OpenMP_CXX-v${OpenMP_CXX_VERSION}")

# MPI
FIND_PACKAGE( MPI REQUIRED )
INCLUDE_DIRECTORIES(SYSTEM ${MPI_CXX_INCLUDE_DIRS})
LIST(APPEND FASTPLI_INFO_LIBRARIES "MPI_CXX-v${MPI_CXX_VERSION}")

# OpenGL
SET(OpenGL_GL_PREFERENCE "GLVND")
FIND_PACKAGE( OpenGL REQUIRED )
FIND_PACKAGE( GLUT REQUIRED )

IF ( NOT OPENGL_GLU_FOUND )
   MESSAGE( SEND_ERROR "OpenGL::GLU not found" )
ENDIF ( NOT OPENGL_GLU_FOUND )

LIST(APPEND FASTPLI_INFO_LIBRARIES "OpenGL-v${OpenGL_VERSION}")
LIST(APPEND FASTPLI_INFO_LIBRARIES "GLUT-v${GLUT_VERSION}")

##########################################################################
# BUILD VERSION FILES
##########################################################################
CONFIGURE_FILE(
   ${CMAKE_SOURCE_DIR}/cmake/init.cmake
   ${CMAKE_CURRENT_BINARY_DIR}/configure_init.cmake
)

ADD_CUSTOM_TARGET(
   init_files ALL
   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/configure_init.cmake
)

##########################################################################
# Sources
##########################################################################
INCLUDE_DIRECTORIES( AFTER src)

SET(PROJECT_LIB_DIR "${PROJECT_SOURCE_DIR}/src/fastpli")
ADD_SUBDIRECTORY("${PROJECT_SOURCE_DIR}/src/model/solver")
ADD_SUBDIRECTORY("${PROJECT_SOURCE_DIR}/src/objects")
ADD_SUBDIRECTORY("${PROJECT_SOURCE_DIR}/src/simulation")
