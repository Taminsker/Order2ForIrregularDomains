
option(EXAMPLE_07 "Build Example_07" ON)

if(EXAMPLE_07)
    message("-- Example_07 ON")
else(EXAMPLE_07)
    message("-- Example_07 OFF")
endif(EXAMPLE_07)

if(EXAMPLE_07)
    add_executable(
        Example_07
        main.cpp
        builder.cpp
        headers.h
        )

    target_link_libraries(Example_07 O2fid)
    target_link_libraries(Example_07 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_07)
