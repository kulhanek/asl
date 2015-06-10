# ==============================================================================
# ASL CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(NETCDF_SRC
        main.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(test-netcdf ${NETCDF_SRC})

TARGET_LINK_LIBRARIES(test-netcdf
                         asl
                         ${NETCDF_CLIB_NAME}
                         ${SCIMAFIC_CLIB_NAME}
                         ${HIPOLY_LIB_NAME}
                         )
