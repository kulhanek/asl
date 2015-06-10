#ifndef NetCDFTrajH
#define NetCDFTrajH
// =============================================================================
// ASL - Amber Support Library
// -----------------------------------------------------------------------------
//    Copyright (C) 2010 Petr Kulhanek (kulhanek@chemi.muni.cz)
//    Copyright (C) 2010 Jakub Stepan (xstepan3@chemi.muni.cz)
//
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with this program; if not, write to the Free Software Foundation, Inc.,
//     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
// =============================================================================

#include <ASLMainHeader.hpp>
#include <SmallString.hpp>
#include <AmberTrajectory.hpp>
#include <netcdf.h>

//---------------------------------------------------------------------------

class CNetCDFTraj;

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------

/// trajectory file IO master class

class ASL_PACKAGE CNetCDFTraj {
public:
    CNetCDFTraj(void);
    ~CNetCDFTraj(void);

// information methods --------------------------------------------------------
    /// is NetCDf file?
    static bool IsNetCDFFile(const CSmallString& name);

// executive methods ----------------------------------------------------------
    /// open trajectory file
    bool Open(const CSmallString& name,ETrajectoryOpenMode mode);

    /// read header
    bool ReadHeader(CAmberTopology* p_top);

    /// write header
    bool WriteHeader(CAmberTopology* p_top);

    /// read trajectory snapshot
    bool ReadSnapshot(CAmberRestart* p_snap);

    /// write trajectory snapshot
    bool WriteSnapshot(CAmberRestart* p_snap);

// section of private data -----------------------------------------------------
private:
    int                     NCID;       // net stream id
    ETrajectoryOpenMode     Mode;

    // header
    CSmallString            Title;
    CSmallString            Application;
    CSmallString            Program;
    CSmallString            Version;
    CSmallString            Conventions;
    CSmallString            ConventionVersion;

    int                     TotalSnapshots;
    int                     ActualAtoms;
    int                     NumOfTopologyAtoms;
    bool                    HasBox;

    int                     SpatialVID;
    int                     SpatialDID;
    int                     Spatial;

    int                     CurrentSnapshot;
    int                     CoordinateVID;
    int                     CoordinateDID;
    float*                  Coordinates;

    int                     CellLengthVID;
    int                     CellLengthDID;
    double                  CellLength[3];

    int                     CellAngleVID;
    int                     CellAngleDID;
    double                  CellAngle[3];

    int                     LabelDID;

    int                     TimeVID;
    int                     TimeDID;
    float                   Time;

    int  GetDimensionInfo(const char* p_attribute, int* p_length);
    bool DefineDimension(const char *name, int length, int *dimidp);
    int  GetVariableID(const char* p_variable);
    bool GetVariableAttribute(int varid,const char* p_attribute,CSmallString& text);
    bool PutAttributeText(int vid, const char *attribute, const char *text);
    bool DefineVariable(const char *name, nc_type xtype, int ndims, int dimids[], int *varidp);

    friend class CAmberTrajectory;
};

//---------------------------------------------------------------------------
#endif
