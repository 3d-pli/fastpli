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

pybind11_add_module(__objects
      bindings/objects_module.cpp
      bindings/fiber_module.cpp
      bindings/cell_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>)

GET_TARGET_PROPERTY(outdir __objects LIBRARY_OUTPUT_DIRECTORY)
IF (outdir)
   SET(outdir "${outdir}/${subdir}")
ELSE()
   SET(outdir "${subdir}")
endif()

SET_TARGET_PROPERTIES(__objects PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${outdir})
