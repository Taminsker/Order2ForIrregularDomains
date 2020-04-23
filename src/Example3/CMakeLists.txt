

option(EXAMPLE_3 "Build Example 3" ON)

if(EXAMPLE_3)
    message("-- Example 3 ON")
else(EXAMPLE_3)
    message("-- Example 3 OFF")
endif(EXAMPLE_3)

if(EXAMPLE_3)
    add_executable(
        Example3
        main.cpp
        )

    target_link_libraries(Example3 O2fid)
    target_link_libraries(Example3 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_3)
