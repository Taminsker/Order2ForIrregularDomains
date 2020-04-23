

option(EXAMPLE_10 "Build Example 10" OFF)

if(EXAMPLE_10)
    message("-- Example 10 ON")
else(EXAMPLE_10)
    message("-- Example 10 OFF")
endif(EXAMPLE_10)

if(EXAMPLE_10)
    add_executable(
        Example10
        main.cpp
        )

    target_link_libraries(Example10 O2fid)
    target_link_libraries(Example10 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_10)
