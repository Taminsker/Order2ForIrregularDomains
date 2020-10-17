
option(EXAMPLE_10 "Build Example_10" OFF)

if(EXAMPLE_10)
    message("-- Example_10 ON")
else(EXAMPLE_10)
    message("-- Example_10 OFF")
endif(EXAMPLE_10)

if(EXAMPLE_10)
    add_executable(
        Example_10
        main.cpp
        builder.cpp
        headers.h
        )

    target_link_libraries(Example_10 O2fid)
    target_link_libraries(Example_10 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_10)
