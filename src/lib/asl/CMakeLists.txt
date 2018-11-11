# ==============================================================================
# ASL CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(ASL_SRC
        ASLMainHeader.cpp

     # types -----------------
        topology/AmberBondType.cpp
        topology/AmberAngleType.cpp
        topology/AmberDihedralType.cpp

     # elements -------------
        topology/AmberAtom.cpp
        topology/AmberBond.cpp
        topology/AmberAngle.cpp
        topology/AmberDihedral.cpp
        topology/AmberResidue.cpp

     # lists ----------------
        topology/AmberAtomList.cpp
        topology/AmberBondList.cpp
        topology/AmberAngleList.cpp
        topology/AmberDihedralList.cpp
        topology/AmberResidueList.cpp
        topology/AmberNonBondedList.cpp

     # special parts --------
        topology/AmberBox.cpp
        topology/AmberCap.cpp

     # topology -------------
        topology/AmberTopology.cpp
        topology/AmberSubTopology.cpp

     # netcdf support
        NetCDFFile.cpp

     # trajectory -----------
        trajectory/AmberTrajectory.cpp
        trajectory/NetCDFTraj.cpp

     # restart --------------
        restart/AmberRestart.cpp
        restart/NetCDFRst.cpp

     # masks ----------------
        mask/AmberMaskAtoms.cpp
        mask/AmberMaskASelection.cpp

        mask/AmberMaskResidues.cpp
        mask/AmberMaskRSelection.cpp

        mask/maskparser/AmberMaskParser.cpp
        mask/maskparser/AmberMask.tab.c
        mask/maskparser/AmberMask.yy.c
        )


# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(asl_static STATIC ${ASL_SRC})
    SET_TARGET_PROPERTIES(asl_static PROPERTIES
                            OUTPUT_NAME asl
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_DEFINITIONS(-DASL_BUILDING_DLL)
    ADD_LIBRARY(asl_shared SHARED ${ASL_SRC})
    TARGET_LINK_LIBRARIES(asl_shared ${SYSTEM_LIBS})
    SET_TARGET_PROPERTIES(asl_shared PROPERTIES
                            OUTPUT_NAME asl
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${ASL_VERSION})
    INSTALL(TARGETS asl_shared
            DESTINATION lib)
ENDIF(LIBS_SHARED)

#-------------------------------------------------------------------------------
