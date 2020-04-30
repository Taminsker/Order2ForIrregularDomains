

option(EXAMPLE_TVD_WENO "Build Example TVD WENO" ON)

if(EXAMPLE_TVD_WENO)
    message("-- Example TVD WENO ON")
else(EXAMPLE_TVD_WENO)
    message("-- Example TVD WENO OFF")
endif(EXAMPLE_TVD_WENO)

if(EXAMPLE_TVD_WENO)
    add_executable(
        ExampleTVD_WENO
        main.cpp
        )

    target_link_libraries(ExampleTVD_WENO O2fid)
    target_link_libraries(ExampleTVD_WENO ${EIGEN3_LIBRARIES})

endif(EXAMPLE_TVD_WENO)
