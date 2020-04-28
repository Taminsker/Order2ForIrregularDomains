

option(EXAMPLE_8C "Build Example 8C" ON)

if(EXAMPLE_8C)
    message("-- Example 8C ON")
else(EXAMPLE_8C)
    message("-- Example 8C OFF")
endif(EXAMPLE_8C)

if(EXAMPLE_8C)
    add_executable(
        Example8C
        main.cpp
        builder.cpp
        headers.h
        )

    target_link_libraries(Example8C O2fid)
    target_link_libraries(Example8C ${EIGEN3_LIBRARIES})

endif(EXAMPLE_8C)
