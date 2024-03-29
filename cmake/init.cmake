set(CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR})
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})
set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(FASTPLI_INFO_COMPILER ${FASTPLI_INFO_COMPILER})
set(FASTPLI_INFO_LIBRARIES ${FASTPLI_INFO_LIBRARIES})

message("configure version files")

if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND git describe --always
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_DESCRIBE_LOG
    OUTPUT_STRIP_TRAILING_WHITESPACE)
else(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  message("did not find git")
  set(GIT_BRANCH "")
  set(GIT_COMMIT_HASH "")
  set(GIT_DESCRIBE_LOG "")
endif(EXISTS "${CMAKE_SOURCE_DIR}/.git")

configure_file(${CMAKE_SOURCE_DIR}/cmake/version.hpp
               ${CMAKE_SOURCE_DIR}/src/include/version.hpp)
configure_file(${CMAKE_SOURCE_DIR}/cmake/__version.py
               ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/__version.py)
configure_file(${CMAKE_SOURCE_DIR}/cmake/setup.py ${CMAKE_SOURCE_DIR}/setup.py)
