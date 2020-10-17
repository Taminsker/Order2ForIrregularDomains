
option(EXAMPLE_12 "Build Example_12" OFF)

if(EXAMPLE_12)
    message("-- Example_12 ON")
else(EXAMPLE_12)
    message("-- Example_12 OFF")
endif(EXAMPLE_12)

if(EXAMPLE_12)
    add_executable(
        Example_12
        main.cpp
        )

    target_link_libraries(Example_12 O2fid)
    target_link_libraries(Example_12 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_12)
