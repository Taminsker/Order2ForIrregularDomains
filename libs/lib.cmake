# --------------------------------------------- #
# ----------------- LIBRARIES ----------------- #
# --------------------------------------------- #

# -------------- BOOST -------------------#
#find_package(Boost REQUIRED COMPONENTS filesystem regex)
#
#if(NOT Boost_FOUND)
#    message(STATUS "find_package : NOT found BOOST library")
#else ()
#    message(STATUS "find_package : found BOOST library")
#    include_directories(${Boost_INCLUDE_DIRS})
#    target_link_libraries(${TARGET} ${Boost_LIBRARIES})
#endif()



# --------------- VTK --------------- #
set(VTK_DIR "${LIB_DIR}/VTK-8.2.0//lib/cmake/vtk-8.2")
set(_vtk_installed_prefix "${CMAKE_SOURCE_DIR}/lib/VTK-8.2.0/lib/cmake/vtk-8.2")

find_package(VTK REQUIRED PATHS ${VTK_DIR})
if (NOT VTK_FOUND)
    message(STATUS "find_package : NOT found VTK library in ${CMAKE_SOURCE_DIR}/lib/VTK-8.2.0/lib/cmake/vtk-8.2 ")
else()
    message(STATUS "find_package : found VTK library in ${CMAKE_SOURCE_DIR}/lib/VTK-8.2.0/lib/cmake/vtk-8.2 ")
    include(${VTK_USE_FILE})
    target_link_libraries(${TARGET} ${VTK_LIBRARIES})
endif()


# --------------- QT --------------- #
#find_package(Qt5 COMPONENTS Widgets REQUIRED)
#if (NOT Qt5_FOUND)
#    message(STATUS "find_package : NOT found Qt5::Widgets library")
#else()
#    message(STATUS "find_package : found Qt5::Widgets library")
#    target_link_libraries(${TARGET} Qt5::Widgets)
#endif()
#ADD_LIBRARY(Qt5::Widgets SHARED IMPORTED ${SRC_FILES})

set(EIGEN3_DIR "${LIB_DIR}/EIGEN/share/eigen3/cmake")

# --------------- EIGEN --------------- #

find_package(Eigen3 REQUIRED PATHS ${EIGEN3_DIR})

if (NOT EIGEN3_FOUND)
    message(STATUS "find_package : NOT found EIGEN library in ${CMAKE_SOURCE_DIR}/lib//EIGEN/share/eigen3/cmake ")
else()
    message(STATUS "find_package : found EIGEN library in ${CMAKE_SOURCE_DIR}/lib//EIGEN/share/eigen3 ")
    include(${EIGEN3_USE_FILE})
    include_directories(${EIGEN3_INCLUDE_DIR})
    include_directories("${EIGEN3_INCLUDE_DIR}/Eigen/unsupported")
    target_link_libraries(${TARGET} ${EIGEN3_LIBRARIES})
endif()


