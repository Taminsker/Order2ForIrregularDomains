

option(BAC_A_SABLE "Build Bac à sable" ON)

if(BAC_A_SABLE)
    message("-- Bac_a_sable ON")
else(BAC_A_SABLE)
    message("-- Bac_a_sable OFF")
endif(BAC_A_SABLE)

if(BAC_A_SABLE)
    add_executable(
        BacASable
        main.cpp
        )

    target_link_libraries(BacASable O2fid)
    target_link_libraries(BacASable ${EIGEN3_LIBRARIES})

endif(BAC_A_SABLE)
