# ==============================================================================
# ASL CMake File
# ==============================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)
#SET(CMAKE_VERBOSE_MAKEFILE ON)
PROJECT(ASL)

# ==============================================================================
# conditional setup ------------------------------------------------------------
# ==============================================================================

SET(LIBS_STATIC OFF CACHE BOOL "Should the static version of hipoly library be built?")
SET(LIBS_SHARED ON CACHE BOOL "Should the dynamic version of hipoly library be built?")
SET(TRY_QT_LIB ON CACHE BOOL "Should the qt lib be used?")

# ==============================================================================
# project setup ----------------------------------------------------------------
# ==============================================================================

# define project name ----------------------------
SET(ASL_VERSION 2.0.x)

# what languages we use --------------------------
ENABLE_LANGUAGE(CXX)

# ==============================================================================
# external dependencies --------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    SET(DEVELOPMENT_ROOT "$ENV{DEVELOPMENT_ROOT}")
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

IF(NOT DEVELOPMENT_ROOT)
    MESSAGE(FATAL_ERROR "DEVELOPMENT_ROOT variable is not set!")
ENDIF(NOT DEVELOPMENT_ROOT)

# ==============================================================================
# external libraries -----------------------------------------------------------
# ==============================================================================

# HIPOLY ====================
SET(HIPOLY_ROOT ${DEVELOPMENT_ROOT}/projects/hipoly/1.0)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/color SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/compress SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/error SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/event SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/file SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/list SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/option SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/string SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/thread SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/time SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/xml SYSTEM)
LINK_DIRECTORIES(${HIPOLY_ROOT}/lib)
SET(HIPOLY_LIB_NAME hipoly)

# SCIMAFIC ====================
SET(SCIMAFIC_ROOT ${DEVELOPMENT_ROOT}/projects/scimafic/1.0)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic SYSTEM)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/maths SYSTEM)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/files SYSTEM)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/chemistry SYSTEM)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/physics SYSTEM)
LINK_DIRECTORIES(${SCIMAFIC_ROOT}/lib)
SET(SCIMAFIC_CLIB_NAME cscimafic)

# NetCDF =====================
SET(NETCDF_ROOT ${DEVELOPMENT_ROOT}/thirdparty/netcdfcore/4.0.1)
INCLUDE_DIRECTORIES(${NETCDF_ROOT}/include SYSTEM)
LINK_DIRECTORIES(${NETCDF_ROOT}/lib)
SET(NETCDF_CLIB_NAME cnetcdf)

SET(SYSTEM_LIBS ${NETCDF_CLIB_NAME}
                ${HIPOLY_LIB_NAME}
                ${SCIMAFIC_CLIB_NAME}
                )

# architecture -----------------------------------
IF(UNIX)
  SET(UNIX 1)
  ADD_DEFINITIONS(-DUNIX)
ENDIF(UNIX)

# policies ---------------------------------------
IF(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

# ==============================================================================
# output directories -----------------------------------------------------------
# ==============================================================================

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
LINK_DIRECTORIES(${LIBRARY_OUTPUT_PATH})

# ==============================================================================
# compiler setups --------------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    # release setup
    SET(CMAKE_C_FLAGS_RELEASE "-O2")
    SET(CMAKE_CXX_FLAGS_RELEASE "-O2")

    # debug setup
    SET(CMAKE_C_FLAGS_DEBUG "-g -O0")
    SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Werror -pedantic-errors -Wundef -Wno-long-long")

    #SET(CMAKE_BUILD_TYPE Release)
    SET(CMAKE_BUILD_TYPE Debug)
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

# ==============================================================================
# project subdirectories  ------------------------------------------------------
# ==============================================================================

ADD_SUBDIRECTORY(src)
