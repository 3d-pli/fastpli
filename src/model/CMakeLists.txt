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
ENDIF ( NOT OPENGL_FOUND OR NOT GLUT_FOUND )
ADD_LIBRARY(${CoreModules} OBJECT ${solverCoreFiles})
SET_PROPERTY(TARGET ${CoreModules} PROPERTY POSITION_INDEPENDENT_CODE TRUE)
IF ( OPENGL_FOUND AND GLUT_FOUND )
   TARGET_COMPILE_DEFINITIONS(${CoreModules} PRIVATE _VIS_LIBRARIES)
ENDIF ( OPENGL_FOUND AND GLUT_FOUND )

# building python modules
pybind11_add_module(__solver
      solver/bindings/solver_module.cpp
      $<TARGET_OBJECTS:${CoreModules}>
      $<TARGET_OBJECTS:objectsCoreModules>
      )
TARGET_LINK_LIBRARIES(__solver PRIVATE ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

GET_TARGET_PROPERTY(outdir __solver LIBRARY_OUTPUT_DIRECTORY)
IF (outdir)
   SET(outdir "${outdir}/${subdir}")
ELSE()
   SET(outdir "${subdir}")
endif()

SET_TARGET_PROPERTIES(__solver PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${outdir})
IF ( OPENGL_FOUND AND GLUT_FOUND )
   TARGET_COMPILE_DEFINITIONS(__solver PRIVATE _VIS_LIBRARIES)
ENDIF ( OPENGL_FOUND AND GLUT_FOUND )
