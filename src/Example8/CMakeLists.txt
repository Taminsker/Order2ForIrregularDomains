

option(EXAMPLE_8 "Build Example 8" ON)

if(EXAMPLE_8)
    message("-- Example 8 ON")
else(EXAMPLE_8)
    message("-- Example 8 OFF")
endif(EXAMPLE_8)

if(EXAMPLE_8)
    add_executable(
        Example8
        main.cpp
        )

    target_link_libraries(Example8 O2fid)
    target_link_libraries(Example8 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_8)
