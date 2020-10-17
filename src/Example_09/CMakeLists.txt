
option(EXAMPLE_09 "Build Example_09" ON)

if(EXAMPLE_09)
    message("-- Example_09 ON")
else(EXAMPLE_09)
    message("-- Example_09 OFF")
endif(EXAMPLE_09)

if(EXAMPLE_09)
    add_executable(
        Example_09
        main.cpp
        builder.cpp
        headers.h
        )

    target_link_libraries(Example_09 O2fid)
    target_link_libraries(Example_09 ${EIGEN3_LIBRARIES})

endif(EXAMPLE_09)
