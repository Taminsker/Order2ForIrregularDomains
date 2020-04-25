

option(EXAMPLE_6 "Build Example 6" ON)

if(EXAMPLE_6)
    message("-- Example 6 ON")
else(EXAMPLE_6)
    message("-- Example 6 OFF")
endif(EXAMPLE_6)

if(EXAMPLE_6)
    add_executable(
        Example6
        main.cpp
        )

    target_link_libraries(Example6 O2fid)
    target_link_libraries(Example6 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_6)
