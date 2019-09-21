CMAKE_MINIMUM_REQUIRED(VERSION 3.9)
SET(subdir simulation)

##########################################################################
# Building shared Library for Python
##########################################################################

# building core objects
SET(CoreModules ${subdir}CoreModules)
FILE(GLOB CoreFiles *.cpp)
ADD_LIBRARY(${CoreModules} OBJECT ${CoreFiles})
SET_PROPERTY(TARGET ${CoreModules} PROPERTY POSITION_INDEPENDENT_CODE TRUE)

# building python modules
pybind11_add_module(__generation
      bindings/generator_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>
      $<TARGET_OBJECTS:objectsCoreModules>
      )
TARGET_LINK_LIBRARIES(__generation PRIVATE ${MPI_LIBRARIES})

pybind11_add_module(__simulation
      bindings/simulator_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>
      $<TARGET_OBJECTS:objectsCoreModules>
      )
TARGET_LINK_LIBRARIES(__simulation PRIVATE ${MPI_LIBRARIES})

SET_TARGET_PROPERTIES(__generation PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}/${subdir}")
SET_TARGET_PROPERTIES(__simulation PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}/${subdir}")

IF ( MPI_FOUND )
   TARGET_COMPILE_DEFINITIONS(${CoreModules} PRIVATE _FASTPLI_MPI_FOUND_)
   TARGET_COMPILE_DEFINITIONS(__generation PRIVATE _FASTPLI_MPI_FOUND_)
   TARGET_COMPILE_DEFINITIONS(__simulation PRIVATE _FASTPLI_MPI_FOUND_)
ENDIF()
