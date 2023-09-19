cmake_minimum_required(VERSION 3.9)
set(subdir model/solver)

# ##############################################################################
# Building shared Library for Python
# ##############################################################################

# building core objects
set(CoreModules solverCoreModules)
file(GLOB solverCoreFiles *.cpp)

if(NOT OPENGL_FOUND OR NOT GLUT_FOUND)
  message(WARNING "OpenGL or GLUT not found. scene.cpp will not be compiled")
  list(REMOVE_ITEM solverCoreFiles ${CMAKE_CURRENT_SOURCE_DIR}/scene.cpp)
endif()

add_library(${CoreModules} OBJECT ${solverCoreFiles})
target_include_directories(${CoreModules} SYSTEM PUBLIC ${PYTHON_INCLUDE_DIRS})
set_property(TARGET ${CoreModules} PROPERTY POSITION_INDEPENDENT_CODE TRUE)

if(OPENGL_FOUND AND GLUT_FOUND)
  target_compile_definitions(${CoreModules} PRIVATE _VIS_LIBRARIES)
endif()

# building python modules
pybind11_add_module(
  __solver bindings/solver_module.cpp $<TARGET_OBJECTS:${CoreModules}>
  $<TARGET_OBJECTS:objectsCoreModules>)

if(OPENGL_FOUND AND GLUT_FOUND)
  target_link_libraries(__solver PRIVATE ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
  target_compile_definitions(__solver PRIVATE _VIS_LIBRARIES)
endif()

set_target_properties(__solver PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                          "${PROJECT_LIB_DIR}/${subdir}")
