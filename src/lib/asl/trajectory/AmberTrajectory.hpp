#ifndef AmberTrajectoryH
#define AmberTrajectoryH
/** \ingroup AmberTrajectory*/
/*! \file AmberTrajectory.hpp */
// =============================================================================
// ASL - Amber Support Library
// -----------------------------------------------------------------------------
//    Copyright (C) 2003,2004,2008 Petr Kulhanek (kulhanek@chemi.muni.cz)
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
#include <stdio.h>
#include <Point.hpp>
#include <SmallString.hpp>

//---------------------------------------------------------------------------

class CAmberTopology;
class CAmberRestart;
class CAmberTrajectory;
class CNetCDFTraj;

//---------------------------------------------------------------------------

/// type of trajectory format
enum ETrajectoryType {
    AMBER_TRAJ_CXYZB, // position + box if IFBOX == 1
    AMBER_TRAJ_VXYZ   // velocities
};

/// open mode for trajectory file
enum ETrajectoryOpenMode {
    AMBER_TRAJ_READ,
    AMBER_TRAJ_WRITE
};

/// trajectory format
enum ETrajectoryFormat {
    AMBER_TRAJ_UNKNOWN,
    AMBER_TRAJ_ASCII,
    AMBER_TRAJ_ASCII_GZIP,
    AMBER_TRAJ_ASCII_BZIP2,
    AMBER_TRAJ_NETCDF
};

//---------------------------------------------------------------------------

/// trajectory file IO master class

class ASL_PACKAGE CAmberTrajectory {
public:
    CAmberTrajectory(void);
    ~CAmberTrajectory(void);

    /// assign topology with trajectory file - previous data are discarded !
    bool AssignTopology(CAmberTopology* p_top);

    /// assign restart with trajectory file - previous data are discarded !
    bool AssignRestart(CAmberRestart* p_rst);

    /// get trajectory info - it will destroy previous contents in class
    bool PrintInfo(const CSmallString& name,ETrajectoryFormat format,
                                ETrajectoryType type,FILE* p_out=NULL);

    /// open trajectory file for reading or writing
    bool OpenTrajectoryFile(const CSmallString& name,
                                         ETrajectoryFormat format,
                                         ETrajectoryType type,
                                         ETrajectoryOpenMode mode);

    /// close trajectory
    bool CloseTrajectoryFile(void);

    /// is it opened?
    bool IsItOpened(void);

    /// assign trajectory to file
    bool AssignTrajectoryToFile(FILE* p_file,
            ETrajectoryFormat format,
            ETrajectoryType type,
            ETrajectoryOpenMode mode);

    /// get title, tittle is read from file in OpenTrajectoryFile or AssignTrajectoryToFile method
    const CSmallString  GetTitle(void);

    /// set title, title is written to file in OpenTrajectoryFile or AssignTrajectoryToFile method
    void  SetTitle(const CSmallString& title);

    /// read snapshot
    /// 0 - OK, 1 - EOF, < 0 - some error
    int ReadSnapshot(void);

    /// read snapshot
    /// 0 - OK, 1 - EOF, < 0 - some error
    int ReadSnapshot(CAmberRestart* p_rst);

    /// write snapshot
    bool WriteSnapshot(void);

    /// write snapshot
    bool WriteSnapshot(CAmberRestart* p_rst);

    /// return snapshot
    CAmberRestart* GetSnapshot(void);

    /// return topology
    CAmberTopology* GetTopology(void);

    /// return number of atoms in trajectory file
    int    GetNumberOfAtoms(void);

    /// return number of snapshots in trajectory (nectcdf) or written to trajectory
    int    GetNumberOfSnapshots(void);

    /// get trajectory format
    ETrajectoryFormat   GetFormat(void);

    /// get trajectory opne mode
    ETrajectoryOpenMode GetOpenMode(void);

// section of private data -----------------------------------------------------
private:
    CAmberTopology*         Topology;
    ETrajectoryType         Type;
    ETrajectoryOpenMode     Mode;
    ETrajectoryFormat       Format;
    FILE*                   TrajectoryFile;
    CNetCDFTraj*            NetCDF;
    bool                    OwnFile;
    bool                    PClose;
    CAmberRestart*          Snapshot;
    char                    Title[81];
    int                     NumOfSnapshots;

    int  ReadSnapshotASCII(CAmberRestart* p_rst);
    bool WriteSnapshotASCII(CAmberRestart* p_rst);
};

//---------------------------------------------------------------------------
#endif
