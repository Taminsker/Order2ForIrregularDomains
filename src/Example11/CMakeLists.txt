

option(EXAMPLE_11 "Build Example 11" OFF)

if(EXAMPLE_11)
    message("-- Example 11 ON")
else(EXAMPLE_11)
    message("-- Example 11 OFF")
endif(EXAMPLE_11)

if(EXAMPLE_11)
    add_executable(
        Example11
        main.cpp
        )

    target_link_libraries(Example11 O2fid)
    target_link_libraries(Example11 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_11)
