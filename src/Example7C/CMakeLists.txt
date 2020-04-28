

option(EXAMPLE_7C "Build Example 7C" ON)

if(EXAMPLE_7C)
    message("-- Example 7C ON")
else(EXAMPLE_7C)
    message("-- Example 7C OFF")
endif(EXAMPLE_7C)

if(EXAMPLE_7C)
    add_executable(
        Example7C
        main.cpp
        builder.cpp
        headers.h
        )

    target_link_libraries(Example7C O2fid)
    target_link_libraries(Example7C ${EIGEN3_LIBRARIES})

endif(EXAMPLE_7C)
