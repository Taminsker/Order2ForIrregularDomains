
option(EXAMPLE_04 "Build Example_04" ON)

if(EXAMPLE_04)
    message("-- Example_04 ON")
else(EXAMPLE_04)
    message("-- Example_04 OFF")
endif(EXAMPLE_04)

if(EXAMPLE_04)
    add_executable(
        Example_04
        main.cpp
        )

    target_link_libraries(Example_04 O2fid)
    target_link_libraries(Example_04 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_04)
