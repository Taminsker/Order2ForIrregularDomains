cmake_minimum_required(VERSION 2.6)

project(Order2ForIrregularDomains)
include(FindPkgConfig)

set(TARGET Order2ForIrregularDomains)

#set(CMAKE_CXX_FLAGS "-Wall -Wall -Wextra -ansi -O -Wwrite-strings -Wuninitialized -Wunreachable-code -Wno-missing-braces -Wno-missing-field-initializers -m64 -O3 -g -std=c++11" )

set(CMAKE_CXX_FLAGS "-Wall -g -std=c++11" )

set(CMAKE_BUILD_TYPE Release)
set(TODAY "13-Jan-2020")

set(CMAKE_BUILD_DIR build)

set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(BUILD_DIR ${CMAKE_SOURCE_DIR}/build)
set(LIB_DIR ${CMAKE_SOURCE_DIR}/libs)


file(GLOB_RECURSE SOURCE_FILES ${SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE HEADER_FILES ${SOURCE_DIR}/*.h)

add_executable(
	${TARGET}
	${SOURCE_FILES}
	${HEADER_FILES}
)

include(${LIB_DIR}/lib.cmake)


#find_package(OpenMP)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(TP PUBLIC OpenMP::OpenMP_CXX)
#endif()

