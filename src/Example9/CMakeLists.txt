

option(EXAMPLE_9 "Build Example 9" ON)

if(EXAMPLE_9)
    message("-- Example 9 ON")
else(EXAMPLE_9)
    message("-- Example 9 OFF")
endif(EXAMPLE_9)

if(EXAMPLE_9)
    add_executable(
        Example9
        main.cpp
        )

    target_link_libraries(Example9 O2fid)
    target_link_libraries(Example9 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_9)
