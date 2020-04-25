

option(EXAMPLE_7 "Build Example 7" ON)

if(EXAMPLE_7)
    message("-- Example 7 ON")
else(EXAMPLE_7)
    message("-- Example 7 OFF")
endif(EXAMPLE_7)

if(EXAMPLE_7)
    add_executable(
        Example7
        main.cpp
        )

    target_link_libraries(Example7 O2fid)
    target_link_libraries(Example7 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_7)
