CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
SET(subdir model)

##########################################################################
# Building shared Library for Python
##########################################################################

# building core objects
SET(CoreModules solverCoreModules)
FILE(GLOB solverCoreFiles solver/*.cpp)
ADD_LIBRARY(${CoreModules} OBJECT ${solverCoreFiles})
SET_PROPERTY(TARGET ${CoreModules} PROPERTY POSITION_INDEPENDENT_CODE TRUE)

# building python modules
pybind11_add_module(_model
      solver/bindings/solver_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>
      $<TARGET_OBJECTS:objectsCoreModules>
      )

GET_TARGET_PROPERTY(outdir _model LIBRARY_OUTPUT_DIRECTORY)
IF (outdir)
   SET(outdir "${outdir}/${subdir}")
ELSE()
   SET(outdir "${subdir}")
endif()

SET_TARGET_PROPERTIES(_model PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${outdir})