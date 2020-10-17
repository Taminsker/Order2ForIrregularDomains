

option(EXAMPLE_03 "Build Example_03" ON)

if(EXAMPLE_03)
    message("-- Example_03 ON")
else(EXAMPLE_03)
    message("-- Example_03 OFF")
endif(EXAMPLE_03)

if(EXAMPLE_03)
    add_executable(
        Example_03
        main.cpp
        )

    target_link_libraries(Example_03 O2fid)
    target_link_libraries(Example_03 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_03)
