
option(EXAMPLE_02 "Build Example_02" ON)

if(EXAMPLE_02)
    message("-- Example_02 ON")
else(EXAMPLE_02)
    message("-- Example_02 OFF")
endif(EXAMPLE_02)

if(EXAMPLE_02)
    add_executable(
        Example_02
        main.cpp
        )

    target_link_libraries(Example_02 O2fid)
    target_link_libraries(Example_02 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_02)
