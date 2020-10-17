

option(EXAMPLE_01 "Build Example_01" ON)

if(EXAMPLE_01)
    message("-- Example_01 ON")
else(EXAMPLE_01)
    message("-- Example_01 OFF")
endif(EXAMPLE_01)

if(EXAMPLE_01)
    add_executable(
        Example_01
        main.cpp
        )

    target_link_libraries(Example_01 O2fid)
    target_link_libraries(Example_01 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_01)
