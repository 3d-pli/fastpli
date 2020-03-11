CMAKE_MINIMUM_REQUIRED(VERSION 3.9)
SET(subdir model/solver)

##########################################################################
# Building shared Library for Python
##########################################################################

# building core objects
SET(CoreModules solverCoreModules)
FILE(GLOB solverCoreFiles *.cpp)
ADD_LIBRARY(${CoreModules} OBJECT ${solverCoreFiles})
SET_PROPERTY(TARGET ${CoreModules} PROPERTY POSITION_INDEPENDENT_CODE TRUE)

# building python modules
pybind11_add_module(__solver
      bindings/solver_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>
      $<TARGET_OBJECTS:objectsCoreModules>
      )

TARGET_LINK_LIBRARIES(__solver PRIVATE ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

SET_TARGET_PROPERTIES(__solver PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                        "${PROJECT_LIB_DIR}/${subdir}")
