

option(EXAMPLE_06 "Build Example_06" ON)

if(EXAMPLE_06)
    message("-- Example_06 ON")
else(EXAMPLE_06)
    message("-- Example_06 OFF")
endif(EXAMPLE_06)

if(EXAMPLE_06)
    add_executable(
        Example_06
        main.cpp
        )

    target_link_libraries(Example_06 O2fid)
    target_link_libraries(Example_06 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_06)
