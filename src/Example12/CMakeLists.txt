

option(EXAMPLE_12 "Build Example 12" OFF)

if(EXAMPLE_12)
    message("-- Example 12 ON")
else(EXAMPLE_12)
    message("-- Example 12 OFF")
endif(EXAMPLE_12)

if(EXAMPLE_12)
    add_executable(
        Example12
        main.cpp
        )

    target_link_libraries(Example12 O2fid)
    target_link_libraries(Example12 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_12)
