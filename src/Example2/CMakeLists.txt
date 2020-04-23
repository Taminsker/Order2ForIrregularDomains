

option(EXAMPLE_2 "Build Example 2" ON)

if(EXAMPLE_2)
    message("-- Example 2 ON")
else(EXAMPLE_2)
    message("-- Example 2 OFF")
endif(EXAMPLE_2)

if(EXAMPLE_2)
    add_executable(
        Example2
        main.cpp
        )

    target_link_libraries(Example2 O2fid)
    target_link_libraries(Example2 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_2)
