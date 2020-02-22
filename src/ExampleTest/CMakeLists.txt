

option(EXAMPLE_TEST "Build Example Test" ON)

if(EXAMPLE_TEST)
    add_executable(
        ExampleTest
        main.cpp
        )

    target_link_libraries(ExampleTest O2fid)
    target_link_libraries(ExampleTest ${EIGEN3_LIBRARIES})

endif()
