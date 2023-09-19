cmake_minimum_required(VERSION 3.9)
set(subdir objects)

# ##############################################################################
# Building shared Library for Python
# ##############################################################################

# building core objects
set(CoreModules ${subdir}CoreModules)
file(GLOB CoreFiles *.cpp)
add_library(${CoreModules} OBJECT ${CoreFiles})
set_property(TARGET ${CoreModules} PROPERTY POSITION_INDEPENDENT_CODE TRUE)
