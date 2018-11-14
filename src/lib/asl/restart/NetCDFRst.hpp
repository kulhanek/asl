#ifndef NetCDFRstH
#define NetCDFRstH
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
#include <AmberTrajectory.hpp>
#include <NetCDFFile.hpp>

//---------------------------------------------------------------------------

/// trajectory file IO master class

class ASL_PACKAGE CNetCDFRst :  public CNetCDFFile {
public:
    CNetCDFRst(void);
    ~CNetCDFRst(void);

// executive methods ----------------------------------------------------------
    /// open trajectory file
    bool Open(const CSmallString& name,char mode);

    /// read header
    bool ReadHeader(CAmberTopology* p_top);

    /// write header
    bool WriteHeader(CAmberTopology* p_top,bool velocities);

    /// read trajectory snapshot
    bool ReadSnapshot(CAmberRestart* p_snap);

    /// write trajectory snapshot
    bool WriteSnapshot(CAmberRestart* p_snap);

// section of private data -----------------------------------------------------
private:
    char                    Mode;  // 'r' - read, 'w' - write

    // header
    CSmallString            Title;
    CSmallString            Application;
    CSmallString            Program;
    CSmallString            Version;
    CSmallString            Conventions;
    CSmallString            ConventionVersion;

    int                     NumOfNetCDFAtoms;
    int                     NumOfTopologyAtoms;
    bool                    HasBox;
    bool                    HasVelocities;

    int                     AtomDID;

    int                     SpatialVID;
    int                     SpatialDID;
    int                     Spatial;

    int                     CoordinateVID;
    double*                 Coordinates;

    int                     VelocityVID;
    double*                 Velocities;

    int                     CellSpatialVID;
    int                     CellSpatialDID;
    int                     CellLengthVID;
    double                  CellLength[3];

    int                     CellAngularVID;
    int                     CellAngularDID;
    int                     CellAngleVID;
    double                  CellAngle[3];

    int                     LabelDID;

    int                     TimeVID;
    int                     TimeDID;
    double                  Time;

    friend class CAmberTrajectory;
};

//---------------------------------------------------------------------------
#endif
