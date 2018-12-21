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

# building python modules
pybind11_add_module(_fiber_cpp
      bindings/fiber_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>)

GET_TARGET_PROPERTY(outdir _fiber_cpp LIBRARY_OUTPUT_DIRECTORY)
IF (outdir)
   SET(outdir "${outdir}/${subdir}")
ELSE()
   SET(outdir "${subdir}")
endif()

SET_TARGET_PROPERTIES(_fiber_cpp PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${outdir})