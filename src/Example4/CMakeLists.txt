

option(EXAMPLE_4 "Build Example 4" ON)

if(EXAMPLE_4)
    message("-- Example 4 ON")
else(EXAMPLE_4)
    message("-- Example 4 OFF")
endif(EXAMPLE_4)

if(EXAMPLE_4)
    add_executable(
        Example4
        main.cpp
        )

    target_link_libraries(Example4 O2fid)
    target_link_libraries(Example4 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_4)
