cmake_minimum_required(VERSION 2.6)

project(Order2ForIrregularDomains)
include(FindPkgConfig)


#set(CMAKE_CXX_FLAGS "-Wall -Wall -Wextra -ansi -O -Wwrite-strings -Wuninitialized -Wunreachable-code -Wno-missing-braces -Wno-missing-field-initializers -m64 -O3 -g -std=c++11" )

set(CMAKE_CXX_FLAGS "-Wall -g -std=c++11" )

set(CMAKE_BUILD_TYPE Release)
set(TODAY "13-Jan-2020")

set(LIB_DIR ${CMAKE_SOURCE_DIR}/libs)

include(${LIB_DIR}/lib.cmake)

add_subdirectory (${CMAKE_SOURCE_DIR}/src/O2FID)
add_subdirectory (${CMAKE_SOURCE_DIR}/src/TestChaleur)
add_subdirectory (${CMAKE_SOURCE_DIR}/src/Example1)
