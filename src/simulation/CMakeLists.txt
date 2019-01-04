CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
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
pybind11_add_module(generation
      bindings/pli_generator_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>)

pybind11_add_module(simulation
      bindings/pli_simulator_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>)

pybind11_add_module(_helper
      bindings/pli_helper_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>)

GET_TARGET_PROPERTY(outdir generation LIBRARY_OUTPUT_DIRECTORY)
IF (outdir)
   SET(outdir "${outdir}/${subdir}")
ELSE()
   SET(outdir "${subdir}")
endif()

SET_TARGET_PROPERTIES(generation PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${outdir})
SET_TARGET_PROPERTIES(simulation PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${outdir})
SET_TARGET_PROPERTIES(_helper PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${outdir})
