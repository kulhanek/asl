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

#include <AmberTrajectory.hpp>
#include <AmberTopology.hpp>
#include <AmberRestart.hpp>
#include "FortranIO.hpp"
#include <string.h>
#include <ErrorSystem.hpp>
#include <errno.h>
#include <FileName.hpp>
#include <NetCDFTraj.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberTrajectory::CAmberTrajectory(void)
{
    Topology = NULL;
    TrajectoryFile = NULL;
    OwnFile = false;
    Snapshot = NULL;
    memset(Title,0,81);
    Format = AMBER_TRAJ_UNKNOWN;
    Mode = AMBER_TRAJ_READ;
    PClose = false;
    NetCDF = NULL;
    NumOfSnapshots = -1;
}

//---------------------------------------------------------------------------

CAmberTrajectory::~CAmberTrajectory(void)
{
    if( (TrajectoryFile != NULL) && (OwnFile == true) ) {
        if( PClose == true ) {
            fclose(TrajectoryFile);
        } else {
            pclose(TrajectoryFile);
        }
    }
    OwnFile = false;
    Snapshot = NULL;
    if( NetCDF != NULL ) delete NetCDF;
    NetCDF = NULL;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTrajectory::AssignTopology(CAmberTopology* p_top)
{
    if( IsItOpened() ) CloseTrajectoryFile();

    TrajectoryFile = NULL;
    Topology = p_top;
    Snapshot = NULL;
    OwnFile = false;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTrajectory::AssignRestart(CAmberRestart* p_rst)
{
    if( IsItOpened() ) CloseTrajectoryFile();

    TrajectoryFile = NULL;
    Snapshot = p_rst;
    if( Topology == NULL ) {
        if( Snapshot != NULL ) Topology = Snapshot->GetTopology();
    }
    OwnFile = false;

    if( (Snapshot != NULL) && (Topology != NULL) ) {
        if( Snapshot->GetNumberOfAtoms() != Topology->AtomList.GetNumberOfAtoms() ) {
            ES_ERROR("number of atoms in topology and restart file differs");
            return(false);
        }
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTrajectory::PrintInfo(const CSmallString& name,
        ETrajectoryFormat format,
        ETrajectoryType type,FILE* p_out)
{
    if( p_out == NULL ) {
        p_out = stdout;
    }

    fprintf(p_out,"\n");

    if( Topology == NULL ) {
        ES_ERROR("topology is not assigned");
        return(false);
    }

    if( OpenTrajectoryFile(name,format,type,AMBER_TRAJ_READ) == false ) {
        ES_ERROR("unable to open trajectory");
        return(false);
    }

    fprintf(p_out," Trajectory name     : %s\n",(const char*)name);
    fprintf(p_out," Trajectory title    : %s\n",(const char*)GetTitle());
    switch(Format) {
    case AMBER_TRAJ_ASCII:
        fprintf(p_out," Format              : ASCII\n");
        break;
    case AMBER_TRAJ_ASCII_GZIP:
        fprintf(p_out," Format              : ASCII (compressed by gzip)\n");
        break;
    case AMBER_TRAJ_ASCII_BZIP2:
        fprintf(p_out," Format              : ASCII (compressed by bzip2)\n");
        break;
    case AMBER_TRAJ_NETCDF:
        fprintf(p_out," Format              : NetCDF\n");
        break;
    case AMBER_TRAJ_UNKNOWN:
    default:
        fprintf(p_out," Format              : unknown\n");
        break;
    }

    int number_of_snapshots = 0;

    if( Format != AMBER_TRAJ_NETCDF ) {
        fprintf(p_out," Number of snapshots : (scanning in progress, please wait)\n");
        while( ReadSnapshot() ) number_of_snapshots++;
    } else {
        number_of_snapshots = NetCDF->TotalSnapshots;
    }

    CloseTrajectoryFile();

    fprintf(p_out," Number of snapshots : %d\n",number_of_snapshots);

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTrajectory::OpenTrajectoryFile(const CSmallString& name,
        ETrajectoryFormat format,
        ETrajectoryType type,
        ETrajectoryOpenMode mode)
{
    if( Topology == NULL ) {
        ES_ERROR("topology is not assigned");
        return(false);
    }

    if( (TrajectoryFile != NULL) || (NetCDF != NULL) ) {
        ES_ERROR("trajectory is already opened");
        return(false);
    }
    NumOfSnapshots = -1;

    // get file name extension
    CFileName file_name(name);

    if( format == AMBER_TRAJ_UNKNOWN ) {
        // try autodetection

        Format = AMBER_TRAJ_ASCII;

        if( file_name.GetFileNameExt() == ".gz" ) {
            Format = AMBER_TRAJ_ASCII_GZIP;
        }
        if( file_name.GetFileNameExt() == ".bz2" ) {
            Format = AMBER_TRAJ_ASCII_BZIP2;
        }

        if( Format == AMBER_TRAJ_ASCII ) {
            if( mode == AMBER_TRAJ_READ ) {
                if( CNetCDFTraj::IsNetCDFFile(name) == true ) {
                    Format = AMBER_TRAJ_NETCDF;
                }
            } else {
                if( file_name.GetFileNameExt() == ".netcdf" ) {
                    Format = AMBER_TRAJ_NETCDF;
                }
            }
        }
    } else {
        Format = format;
    }

    FILE* p_trajfile;
    PClose = false;

    switch(Format) {
    case AMBER_TRAJ_ASCII:
        if( mode == AMBER_TRAJ_READ ) {
            p_trajfile = fopen(name,"rb");
        } else {
            p_trajfile = fopen(name,"wb");
        }
        break;
    case AMBER_TRAJ_ASCII_GZIP: {
        if( mode == AMBER_TRAJ_READ ) {
            CSmallString command;
            command << "gzip -d -c '" << name << "'";
            p_trajfile = popen(command,"r");
            PClose = true;
        } else {
            CSmallString command;
            command << "gzip -c '" << name << "'";
            p_trajfile = popen(command,"w");
            PClose = true;
        }
    }
    break;
    case AMBER_TRAJ_ASCII_BZIP2: {
        if( mode == AMBER_TRAJ_READ ) {
            CSmallString command;
            command << "bunzip2 -c '" << name << "'";
            p_trajfile = popen(command,"r");
            PClose = true;
        } else {
            CSmallString command;
            command << "bzip2 -c '" << name << "'";
            p_trajfile = popen(command,"w");
            PClose = true;
        }
    }
    break;
    case AMBER_TRAJ_NETCDF:
        NetCDF = new CNetCDFTraj();
        if( NetCDF->Open(name,mode) == false ){
            ES_TRACE_ERROR("unable to open NetCDF");
            return(false);
        }
        if( mode == AMBER_TRAJ_READ ) {
            if( NetCDF->ReadHeader(Topology) == false ){
                ES_TRACE_ERROR("unable to read header");
                return(false);
            }
            NumOfSnapshots = NetCDF->TotalSnapshots;
            Mode = AMBER_TRAJ_READ;
        } else {
            if( NetCDF->WriteHeader(Topology) == false ) {
                ES_TRACE_ERROR("unable to write header");
                return(false);
            }
            NumOfSnapshots = 0;
            Mode = AMBER_TRAJ_WRITE;
        }
        return(true);
    default:
        ES_ERROR("not implemented or supported format");
        return(false);
    }

    if( p_trajfile == NULL ) {
        CSmallString error;
        error << "unable to open file '" << name << "' (" << strerror(errno) << ")";
        ES_ERROR(error);
        return(false);
    }

    OwnFile = true;

    if( AssignTrajectoryToFile(p_trajfile,Format,type,mode) == false ) {
        OwnFile = false;
        fclose(p_trajfile);
        ES_TRACE_ERROR("unable to assign trajectory file");
        return(false);
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberTrajectory::CloseTrajectoryFile(void)
{
    if( NetCDF != NULL ) {
        delete NetCDF;
        NetCDF = NULL;
    }

    if( (TrajectoryFile != NULL) && (OwnFile == true) ) {
        if( PClose == true ) {
            fclose(TrajectoryFile);
        } else {
            pclose(TrajectoryFile);
        }
    }
    TrajectoryFile = NULL;
    NumOfSnapshots = -1;
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberTrajectory::IsItOpened(void)
{
    return( (TrajectoryFile != NULL) || (NetCDF != NULL) );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTrajectory::AssignTrajectoryToFile(FILE* p_file,
        ETrajectoryFormat format,
        ETrajectoryType type,
        ETrajectoryOpenMode mode)
{
    NumOfSnapshots = -1;

    if( NetCDF != NULL ) {
        delete NetCDF;
        NetCDF = NULL;
    }

    if( Topology == NULL ) return(false);

    TrajectoryFile = p_file;
    if( TrajectoryFile == NULL ) return(false);

    if( (format == AMBER_TRAJ_UNKNOWN) ||
            (format == AMBER_TRAJ_NETCDF) ) {
        ES_ERROR("UNKNOWN or NETCDF are not supported formats");
        return(false);
    }

    Format = format;
    Mode = mode;
    Type = type;

    if( Mode == AMBER_TRAJ_READ ) {
        // read title
        CFortranIO fortranio(TrajectoryFile);
        fortranio.SetFormat("1A80");
        if( fortranio.ReadString(Title) == false ) return(false);
    }

    if( Mode == AMBER_TRAJ_WRITE ) {
        // read title
        CFortranIO fortranio(TrajectoryFile);
        fortranio.SetFormat("1A80");
        if( fortranio.WriteString(Title) == false ) return(false);
        fortranio.WriteEndOfSection();
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

const CSmallString CAmberTrajectory::GetTitle(void)
{
    if( NetCDF != NULL ) return(NetCDF->Title);
    return(Title);
}

//---------------------------------------------------------------------------

void CAmberTrajectory::SetTitle(const CSmallString& title)
{
    strncpy(Title,title,80);
}

//------------------------------------------------------------------------------

ETrajectoryFormat CAmberTrajectory::GetFormat(void)
{
    return(Format);
}

//------------------------------------------------------------------------------

ETrajectoryOpenMode CAmberTrajectory::GetOpenMode(void)
{
    return(Mode);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTrajectory::ReadSnapshot(void)
{
    if( Snapshot == NULL ) {
        ES_ERROR("Snapshot is NULL");
        return(false);
    }

    if( NetCDF != NULL ) {
        return(NetCDF->ReadSnapshot(Snapshot));
    } else {
        return( ReadSnapshotASCII(Snapshot) );
    }
}

//---------------------------------------------------------------------------

bool CAmberTrajectory::ReadSnapshot(CAmberRestart* p_rst)
{
    if( Topology == NULL ) {
        ES_ERROR("Topology is NULL");
        return(false);
    }
    if( p_rst == NULL ) {
        ES_ERROR("p_rst is NULL");
        return(false);
    }
    if( p_rst->GetTopology()->AtomList.GetNumberOfAtoms() != Topology->AtomList.GetNumberOfAtoms() ) {
        ES_ERROR("incompatible number of atoms in restart and topology");
        return(false);
    }

    if( p_rst->GetNumberOfAtoms() == 0 ) {
        if( p_rst->Create() == false ) {
            ES_ERROR("unable to initialize snapshot");
            return(false);
        }
    }

    if( NetCDF != NULL ) {
        return(NetCDF->ReadSnapshot(p_rst));
    } else {
        return( ReadSnapshotASCII(p_rst) );
    }
}

//---------------------------------------------------------------------------

bool CAmberTrajectory::WriteSnapshot(void)
{
    if( Snapshot == NULL ) {
        ES_ERROR("Snapshot is NULL");
        return(false);
    }
    bool result;
    if( NetCDF != NULL ) {
        result = NetCDF->WriteSnapshot(Snapshot);
    } else {
        result = WriteSnapshotASCII(Snapshot);
    }
    if( result == true ){
        NumOfSnapshots++;
    }
    return(result);
}

//---------------------------------------------------------------------------

bool CAmberTrajectory::WriteSnapshot(CAmberRestart* p_rst)
{
    if( Topology == NULL ) {
        ES_ERROR("Topology is NULL");
        return(false);
    }
    if( p_rst == NULL ) {
        ES_ERROR("p_rst is NULL");
        return(false);
    }
    if( p_rst->GetTopology()->AtomList.GetNumberOfAtoms() != Topology->AtomList.GetNumberOfAtoms() ) {
        ES_ERROR("incompatible number of atoms in restart and topology");
        return(false);
    }
    if( p_rst->GetNumberOfAtoms() == 0 ) {
        ES_ERROR("restart does not contain any data");
        return(false);
    }
    bool result;
    if( NetCDF != NULL ) {
        result = NetCDF->WriteSnapshot(p_rst);
    } else {
        result = WriteSnapshotASCII(p_rst);
    }
    if( result == true ){
        NumOfSnapshots++;
    }
    return(result);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberRestart* CAmberTrajectory::GetSnapshot(void)
{
    return(Snapshot);
}

//------------------------------------------------------------------------------

CAmberTopology* CAmberTrajectory::GetTopology(void)
{
    return(Topology);
}

//------------------------------------------------------------------------------

int CAmberTrajectory::GetNumberOfAtoms(void)
{
    if( Topology == NULL ) return(-1);
    return( Topology->AtomList.GetNumberOfAtoms() );
}

//------------------------------------------------------------------------------

int CAmberTrajectory::GetNumberOfSnapshots(void)
{
    if( Topology == NULL ) return(-1);
    return( NumOfSnapshots );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTrajectory::ReadSnapshotASCII(CAmberRestart* p_rst)
{
    if( Topology == NULL ) return(false);
    if( p_rst == NULL ) return(false);
    if( TrajectoryFile == NULL ) {
        ES_ERROR("TrajectoryFile is NULL");
        return(false);
    }

    CFortranIO fortranio(TrajectoryFile,true);

    fortranio.SetFormat("10E8.3");

    CPoint* p_data = p_rst->Positions;
    if( p_data == NULL ) return(false);
    bool result = true;

    for(int i=0; i < p_rst->GetNumberOfAtoms(); i++) {
        result &= fortranio.ReadReal(p_data->x);
        result &= fortranio.ReadReal(p_data->y);
        result &= fortranio.ReadReal(p_data->z);
        if( result == false ) return(false);
        p_data++;
    }
    if( Type == AMBER_TRAJ_VXYZ ) return(true);
    if( Topology->BoxInfo.GetType() == AMBER_BOX_NONE ) return(true);

    fortranio.SetFormat("3E8.3");

    result &= fortranio.ReadReal(p_rst->Box.x);
    result &= fortranio.ReadReal(p_rst->Box.y);
    result &= fortranio.ReadReal(p_rst->Box.z);

    return(result);
}

//------------------------------------------------------------------------------

bool CAmberTrajectory::WriteSnapshotASCII(CAmberRestart* p_rst)
{
    if( Topology == NULL ) return(false);
    if( p_rst == NULL ) return(false);
    if( TrajectoryFile == NULL ) {
        ES_ERROR("TrajectoryFile is NULL");
        return(false);
    }

    CFortranIO fortranio(TrajectoryFile,true);

    fortranio.SetFormat("10F8.3");

    CPoint* p_data = p_rst->Positions;
    bool result = true;

    for(int i=0; i < p_rst->GetNumberOfAtoms(); i++) {
        result &= fortranio.WriteReal(p_data->x);
        result &= fortranio.WriteReal(p_data->y);
        result &= fortranio.WriteReal(p_data->z);
        if( result == false ) return(false);
        p_data++;
    }

    fortranio.WriteEndOfSection();

    if( Type == AMBER_TRAJ_VXYZ ) return(true);
    if( Topology->BoxInfo.GetType() == AMBER_BOX_NONE ) return(true);

    fortranio.SetFormat("3F8.3");

    result &= fortranio.WriteReal(p_rst->Box.x);
    result &= fortranio.WriteReal(p_rst->Box.y);
    result &= fortranio.WriteReal(p_rst->Box.z);

    fortranio.WriteEndOfSection();

    return(result);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



