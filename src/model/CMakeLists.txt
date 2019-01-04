CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
SET(subdir model)

##########################################################################
# Building shared Library for Python
##########################################################################

# building core objects
SET(CoreModules solverCoreModules)
FILE(GLOB solverCoreFiles solver/*.cpp)
IF ( NOT OPENGL_FOUND OR NOT GLUT_FOUND )
   LIST(REMOVE_ITEM solverCoreFiles ${CMAKE_CURRENT_SOURCE_DIR}/solver/scene.cpp)
ENDIF ( OPENGL_FOUND AND GLUT_FOUND )
ADD_LIBRARY(${CoreModules} OBJECT ${solverCoreFiles})
SET_PROPERTY(TARGET ${CoreModules} PROPERTY POSITION_INDEPENDENT_CODE TRUE)

# building python modules
pybind11_add_module(_solver_cpp
      solver/bindings/solver_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>
      $<TARGET_OBJECTS:objectsCoreModules>
      )
TARGET_LINK_LIBRARIES(_solver_cpp PRIVATE ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

GET_TARGET_PROPERTY(outdir _solver_cpp LIBRARY_OUTPUT_DIRECTORY)
IF (outdir)
   SET(outdir "${outdir}/${subdir}")
ELSE()
   SET(outdir "${subdir}")
endif()

SET_TARGET_PROPERTIES(_solver_cpp PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${outdir})

IF ( OPENGL_FOUND AND GLUT_FOUND )
      TARGET_COMPILE_DEFINITIONS(${CoreModules} PRIVATE _VIS_LIBRARIES)
      TARGET_COMPILE_DEFINITIONS(_solver_cpp PRIVATE _VIS_LIBRARIES)
# ELSE ( OPENGL_FOUND AND GLUT_FOUND )
#       LIST(REMOVE_ITEM solverCoreFiles ${CMAKE_CURRENT_SOURCE_DIR}/solver/scene.cpp)
ENDIF ( OPENGL_FOUND AND GLUT_FOUND )