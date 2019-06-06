CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
SET(subdir objects)

##########################################################################
# Building shared Library for Python
##########################################################################

# building core objects
SET(CoreModules ${subdir}CoreModules)
FILE(GLOB CoreFiles *.cpp)
ADD_LIBRARY(${CoreModules} OBJECT ${CoreFiles})
SET_PROPERTY(TARGET ${CoreModules} PROPERTY POSITION_INDEPENDENT_CODE TRUE)

pybind11_add_module(__fiber
      bindings/fiber_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>)

pybind11_add_module(__cell
      bindings/cell_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>)

SET_TARGET_PROPERTIES(__fiber PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fastpli/${subdir}")
SET_TARGET_PROPERTIES(__cell PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fastpli/${subdir}")
