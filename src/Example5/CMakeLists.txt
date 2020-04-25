

option(EXAMPLE_5 "Build Example 5" ON)

if(EXAMPLE_5)
    message("-- Example 5 ON")
else(EXAMPLE_5)
    message("-- Example 5 OFF")
endif(EXAMPLE_5)

if(EXAMPLE_5)
    add_executable(
        Example5
        main.cpp
        )

    target_link_libraries(Example5 O2fid)
    target_link_libraries(Example5 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_5)
