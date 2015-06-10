# ==============================================================================
# ASL CMake File
# ==============================================================================

# this is required for asl library ------------------------------------------
INCLUDE_DIRECTORIES(lib/asl SYSTEM)
INCLUDE_DIRECTORIES(lib/asl/topology SYSTEM)
INCLUDE_DIRECTORIES(lib/asl/restart SYSTEM)
INCLUDE_DIRECTORIES(lib/asl/trajectory SYSTEM)
INCLUDE_DIRECTORIES(lib/asl/mask SYSTEM)

# include subdirectories -------------------------------------------------------
ADD_SUBDIRECTORY(lib)
# ADD_SUBDIRECTORY(test)
