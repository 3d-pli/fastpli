CMAKE_MINIMUM_REQUIRED(VERSION 3.9)
SET(subdir simulation)

##########################################################################
# Building shared Library for Python
##########################################################################

# building core objects
SET(CoreModules ${subdir}CoreModules)
ADD_LIBRARY(${CoreModules} OBJECT cell_population.cpp
                                    fiber_bundle.cpp
                                    my_mpi.cpp)
SET_PROPERTY(TARGET ${CoreModules} PROPERTY POSITION_INDEPENDENT_CODE TRUE)
TARGET_INCLUDE_DIRECTORIES(${CoreModules} SYSTEM PUBLIC ${MPI_CXX_INCLUDE_DIRS})

IF ( THESIS )
   TARGET_COMPILE_DEFINITIONS(${CoreModules} PRIVATE _THESIS)
ENDIF ( )

# building python modules
pybind11_add_module(__generation
      bindings/generator_module.cpp
      generator.cpp
      $<TARGET_OBJECTS:${CoreModules}>
      $<TARGET_OBJECTS:objectsCoreModules>
      )
TARGET_LINK_LIBRARIES(__generation PRIVATE ${MPI_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(__generation SYSTEM PUBLIC ${MPI_CXX_INCLUDE_DIRS})

pybind11_add_module(__simulation
      bindings/simulator_module.cpp
      simulator.cpp
      $<TARGET_OBJECTS:${CoreModules}>
      $<TARGET_OBJECTS:objectsCoreModules>
      )
TARGET_LINK_LIBRARIES(__simulation PRIVATE ${MPI_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(__simulation SYSTEM PUBLIC ${MPI_CXX_INCLUDE_DIRS})

# copy to fastpli build python path
SET_TARGET_PROPERTIES(__generation PROPERTIES LIBRARY_OUTPUT_DIRECTORY 
                        "${PROJECT_LIB_DIR}/${subdir}")
SET_TARGET_PROPERTIES(__simulation PROPERTIES LIBRARY_OUTPUT_DIRECTORY 
                        "${PROJECT_LIB_DIR}/${subdir}")

IF ( THESIS )
      TARGET_COMPILE_DEFINITIONS(__simulation PRIVATE _THESIS)
ENDIF ( )
