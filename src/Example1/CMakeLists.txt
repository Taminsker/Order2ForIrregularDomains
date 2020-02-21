

option(EXAMPLE_1 "Build Example 1" ON)

if(EXAMPLE_1)
    add_executable(
        Example1
        main.cpp
        )

    target_link_libraries(Example1 O2fid)
    target_link_libraries(Example1 ${EIGEN3_LIBRARIES})

endif()
