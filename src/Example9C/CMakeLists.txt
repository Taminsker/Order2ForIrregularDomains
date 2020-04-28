

option(EXAMPLE_9C "Build Example 9C" ON)

if(EXAMPLE_9C)
    message("-- Example 9C ON")
else(EXAMPLE_9C)
    message("-- Example 9C OFF")
endif(EXAMPLE_9C)

if(EXAMPLE_9C)
    add_executable(
        Example9C
        main.cpp
        builder.cpp
        headers.h
        )

    target_link_libraries(Example9C O2fid)
    target_link_libraries(Example9C ${EIGEN3_LIBRARIES})

endif(EXAMPLE_9C)
