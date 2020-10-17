
option(EXAMPLE_08 "Build Example_08" ON)

if(EXAMPLE_08)
    message("-- Example_08 ON")
else(EXAMPLE_08)
    message("-- Example_08 OFF")
endif(EXAMPLE_08)

if(EXAMPLE_08)
    add_executable(
        Example_08
        main.cpp
        builder.cpp
        headers.h
        )

    target_link_libraries(Example_08 O2fid)
    target_link_libraries(Example_08 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_08)
