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
