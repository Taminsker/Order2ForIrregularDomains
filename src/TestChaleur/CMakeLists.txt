
option(TEST_CHALEUR "Build TestChaleur" OFF)

if(TEST_CHALEUR)
    add_executable(
        TestChaleur
        main.cpp
        Chaleur1D.cpp
        Chaleur1D.h
        )

    target_link_libraries(TestChaleur O2fid)
    target_link_libraries(TestChaleur ${EIGEN3_LIBRARIES})
endif()
