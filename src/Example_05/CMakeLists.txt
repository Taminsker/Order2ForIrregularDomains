
option(EXAMPLE_05 "Build Example_05" ON)

if(EXAMPLE_05)
    message("-- Example_05 ON")
else(EXAMPLE_05)
    message("-- Example_05 OFF")
endif(EXAMPLE_05)

if(EXAMPLE_05)
    add_executable(
        Example_05
        main.cpp
        )

    target_link_libraries(Example_05 O2fid)
    target_link_libraries(Example_05 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_05)
