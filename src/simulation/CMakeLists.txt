cmake_minimum_required(VERSION 3.9)
set(subdir simulation)

# ##############################################################################
# Building shared Library for Python
# ##############################################################################

# building core objects
set(CoreModules ${subdir}CoreModules)
add_library(${CoreModules} OBJECT cell_population.cpp fiber_bundle.cpp
                                  my_mpi.cpp)
set_property(TARGET ${CoreModules} PROPERTY POSITION_INDEPENDENT_CODE TRUE)
target_include_directories(${CoreModules} SYSTEM PUBLIC ${MPI_CXX_INCLUDE_DIRS})

if(THESIS)
  target_compile_definitions(${CoreModules} PRIVATE _THESIS)
endif()

# building python modules
pybind11_add_module(
  __generation bindings/generator_module.cpp generator.cpp
  $<TARGET_OBJECTS:${CoreModules}> $<TARGET_OBJECTS:objectsCoreModules>)
target_link_libraries(__generation PRIVATE ${MPI_LIBRARIES})
target_include_directories(__generation SYSTEM PUBLIC ${MPI_CXX_INCLUDE_DIRS})

pybind11_add_module(
  __simulation bindings/simulator_module.cpp simulator.cpp
  $<TARGET_OBJECTS:${CoreModules}> $<TARGET_OBJECTS:objectsCoreModules>)
target_link_libraries(__simulation PRIVATE ${MPI_LIBRARIES})
target_include_directories(__simulation SYSTEM PUBLIC ${MPI_CXX_INCLUDE_DIRS})

# copy to fastpli build python path
set_target_properties(__generation PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                              "${PROJECT_LIB_DIR}/${subdir}")
set_target_properties(__simulation PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                              "${PROJECT_LIB_DIR}/${subdir}")

if(THESIS)
  target_compile_definitions(__simulation PRIVATE _THESIS)
endif()
