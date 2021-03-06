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
// this file is adopted from ptraj
// used files: netcdf_ptraj.c netcdf_ptraj.h and trajectory.h
//------------------------------------------------------------------------------

#include <NetCDFTraj.hpp>
#include <ErrorSystem.hpp>
#include <AmberRestart.hpp>
#include <AmberTopology.hpp>
#include <string.h>

#define AMBER_NETCDF_FRAME "frame"
#define AMBER_NETCDF_SPATIAL "spatial"
#define AMBER_NETCDF_ATOM "atom"
#define AMBER_NETCDF_CELL_SPATIAL "cell_spatial"
#define AMBER_NETCDF_CELL_ANGULAR "cell_angular"
#define AMBER_NETCDF_COORDS "coordinates"
#define AMBER_NETCDF_TIME "time"
#define AMBER_NETCDF_LABEL "label"
#define AMBER_NETCDF_LABELLEN 5

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CNetCDFTraj::CNetCDFTraj(void)
{
    Mode = AMBER_TRAJ_READ;

    SpatialVID = -1;
    SpatialDID = -1;

    TotalSnapshots = 0;
    Spatial = 0;
    ActualAtoms = 0;

    CurrentSnapshot = 0;
    CoordinateVID = -1;
    CoordinateDID = -1;
    Coordinates = NULL;

    CellSpatialVID = -1;
    CellSpatialDID = -1;
    CellAngularVID = -1;
    CellAngularDID = -1;

    CellLengthVID = -1;
    CellAngleVID = -1;

    TimeVID = -1;
    TimeDID = -1;
    Time = 0;

    HasBox = false;
    NumOfTopologyAtoms = 0;
}

//---------------------------------------------------------------------------

CNetCDFTraj::~CNetCDFTraj(void)
{
    if( Coordinates != NULL ) delete[] Coordinates;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CNetCDFTraj::Open(const CSmallString& name,ETrajectoryOpenMode mode)
{
    if( NCID >= 0 ) {
        ES_ERROR("file is already opened");
        return(false);
    }
    Mode = mode;

    if( Mode == AMBER_TRAJ_READ ) {
        return(CNetCDFFile::Open(name,'r'));
    }

    if( Mode == AMBER_TRAJ_WRITE ) {
        return(CNetCDFFile::Open(name,'w'));
    }

    ES_ERROR("unsupported mode");
    return(false);
}

//------------------------------------------------------------------------------

bool CNetCDFTraj::ReadHeader(CAmberTopology* p_top)
{
    if( p_top == NULL ){
        INVALID_ARGUMENT("p_top == NULL");
    }

    if( NCID < 0 ) {
        ES_ERROR("file is not opened");
        return(false);
    }

    HasBox = p_top->BoxInfo.GetType() != AMBER_BOX_NONE;
    NumOfTopologyAtoms = p_top->AtomList.GetNumberOfAtoms();

    bool result = true;

    result &= GetVariableAttribute(NC_GLOBAL,"title",Title);
    result &= GetVariableAttribute(NC_GLOBAL,"application",Application);
    result &= GetVariableAttribute(NC_GLOBAL,"program",Program);
    result &= GetVariableAttribute(NC_GLOBAL,"programVersion",Version);
    result &= GetVariableAttribute(NC_GLOBAL,"Conventions",Conventions);
    result &= GetVariableAttribute(NC_GLOBAL,"ConventionVersion",ConventionVersion);

    if( Conventions != "AMBER" ) {
        CSmallString error;
        error << "illegal conventions '" << Conventions << "', expecting 'AMBER'";
        ES_ERROR(error);
        return(false);
    }
    if( ConventionVersion != "1.0" ) {
        CSmallString error;
        error << "illegal convention version '" << ConventionVersion << "', expecting '1.0'";
        ES_ERROR(error);
        return(false);
    }

    GetDimensionInfo("frame", &TotalSnapshots);
    GetDimensionInfo("spatial", &Spatial);
    GetDimensionInfo("atom", &ActualAtoms);

    if( Spatial != 3 ) {
        CSmallString error;
        error << "three dim expected but '" << Spatial << "' provided";
        ES_ERROR(error);
        return(false);
    }

    if( ActualAtoms != NumOfTopologyAtoms ) {
        CSmallString error;
        error << "number of atoms in the topology '" << ActualAtoms << "' is different than in topology '" << NumOfTopologyAtoms << "'";
        ES_ERROR(error);
        return(false);
    }

    // sanity check - spatial variable ---------------------

    SpatialVID = GetVariableID("spatial");
    if( SpatialVID < 0 ) {
        return(false);
    }

    size_t  start[3],count[3];
    char    xyz[3];
    int     err;

    start[0] = 0;
    count[0] = Spatial;
    xyz[0] = (char) 0;

    err = nc_get_vara_text(NCID,SpatialVID,start,count,xyz);
    if( err != NC_NOERR ) {
        CSmallString error;
        error << "unable to get spatial names (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    if( (xyz[0] != 'x') || (xyz[1] != 'y') || (xyz[2] != 'z') ) {
        CSmallString error;
        error << "incorrect spatial labels (" << xyz[0] << "," << xyz[1] << "," << xyz[2] << ")";
        ES_ERROR(error);
        return(false);
    }

    CSmallString unit;

    // sanity check - time ---------------------------------
    // time is optional, it is not in trajectories generated by clustering ...
    // RT#694968
    TimeVID = GetVariableID("time",false);
    if( TimeVID >= 0 ) {
        if( GetVariableAttribute(TimeVID,"units",unit) == false ) {
            return(false);
        }
        if( unit != "picosecond" ) {
            CSmallString error;
            error << "incorrect unit for time (" << unit << "), requested picosecond";
            ES_ERROR(error);
            return(false);
        }
    }

    // sanity check - coordinates --------------------------
    CoordinateVID = GetVariableID("coordinates");
    if( CoordinateVID < 0 ) {
        return(false);
    }
    if( GetVariableAttribute(CoordinateVID,"units",unit) == false ) {
        return(false);
    }
    if( unit != "angstrom" ) {
        CSmallString error;
        error << "incorrect unit for coordinates (" << unit << "), requested angstrom";
        ES_ERROR(error);
        return(false);
    }

    CurrentSnapshot = 0;
    if( Coordinates != NULL ) {
        delete Coordinates;
    }
    Coordinates = new float[ActualAtoms*3];

    // sanity check - box ----------------------------------
    if( HasBox ){
        // optional
        CellLengthVID = GetVariableID("cell_lengths");
        CellAngleVID = GetVariableID("cell_angles");
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CNetCDFTraj::WriteHeader(CAmberTopology* p_top)
{
    if( p_top == NULL ){
        INVALID_ARGUMENT("p_top == NULL");
    }

    if( NCID < 0 ) {
        ES_ERROR("file is not opened");
        return(false);
    }

    int dimensionID[NC_MAX_VAR_DIMS];

    NumOfTopologyAtoms = p_top->AtomList.GetNumberOfAtoms();
    ActualAtoms = NumOfTopologyAtoms;

    CurrentSnapshot = 0;
    TotalSnapshots = 0;
    if( Coordinates != NULL ) {
        delete Coordinates;
    }
    Coordinates = new float[ActualAtoms*3];

    // global dimmensions
    DefineDimension(AMBER_NETCDF_FRAME, NC_UNLIMITED, &TimeDID);
    DefineDimension(AMBER_NETCDF_SPATIAL, 3, &SpatialDID);
    DefineDimension(AMBER_NETCDF_ATOM, ActualAtoms, &CoordinateDID);
    DefineDimension(AMBER_NETCDF_LABEL, AMBER_NETCDF_LABELLEN, &LabelDID);
    DefineDimension(AMBER_NETCDF_CELL_SPATIAL, 3, &CellSpatialDID);
    DefineDimension(AMBER_NETCDF_CELL_ANGULAR, 3, &CellAngularDID);

    // put global attributes
    Conventions =  "AMBER";
    ConventionVersion = "1.0";

    PutAttributeText(NC_GLOBAL, "title", Title);
    PutAttributeText(NC_GLOBAL, "application", Application);
    PutAttributeText(NC_GLOBAL, "program", Program);
    PutAttributeText(NC_GLOBAL, "programVersion", Version);
    PutAttributeText(NC_GLOBAL, "Conventions", Conventions);
    PutAttributeText(NC_GLOBAL, "ConventionVersion", ConventionVersion);

    // handle definition of non-optional variables
    dimensionID[0] = SpatialDID;
    DefineVariable(AMBER_NETCDF_SPATIAL, NC_CHAR, 1, dimensionID, &SpatialVID);

    dimensionID[0] = TimeDID;
    DefineVariable(AMBER_NETCDF_TIME, NC_FLOAT, 1, dimensionID, &TimeVID);
    PutAttributeText(TimeVID, "units", "picosecond");

    dimensionID[0] = TimeDID;
    dimensionID[1] = CoordinateDID;
    dimensionID[2] = SpatialDID;
    DefineVariable(AMBER_NETCDF_COORDS, NC_FLOAT, 3, dimensionID, &CoordinateVID);
    PutAttributeText(CoordinateVID, "units", "angstrom");

    dimensionID[0] = CellSpatialDID;
    DefineVariable(AMBER_NETCDF_CELL_SPATIAL, NC_CHAR, 1, dimensionID, &CellSpatialVID);

    dimensionID[0] = CellAngularDID;
    dimensionID[1] = LabelDID;
    DefineVariable(AMBER_NETCDF_CELL_ANGULAR, NC_CHAR, 2, dimensionID, &CellAngularVID);

    // set up box coords
    HasBox = p_top->BoxInfo.GetType() != AMBER_BOX_NONE;

    if( HasBox ) {
        dimensionID[0] = TimeDID;
        dimensionID[1] = CellSpatialDID;
        DefineVariable("cell_lengths", NC_DOUBLE, 2, dimensionID, &CellLengthVID);
        PutAttributeText(CellLengthVID, "units", "angstrom");

        dimensionID[1] = CellAngularDID;
        DefineVariable("cell_angles", NC_DOUBLE, 2, dimensionID, &CellAngleVID);
        PutAttributeText(CellAngleVID, "units", "degree");
    }

    int err,oldMode;

    // set fill mode
    err = nc_set_fill(NCID, NC_NOFILL, &oldMode);
    if (err != NC_NOERR) {
        CSmallString error;
        error << "unable to set fill value (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    // end definition
    err = nc_enddef(NCID);
    if (err != NC_NOERR) {
        CSmallString error;
        error << "unable to end definitions (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    size_t start[3];
    size_t count[3];
    char xyz[3];
    char abc[15] = { 'a', 'l', 'p', 'h', 'a',
                     'b', 'e', 't', 'a', ' ',
                     'g', 'a', 'm', 'm', 'a' };

    // setup labels
    start[0] = 0;
    count[0] = 3;
    xyz[0] = 'x';
    xyz[1] = 'y';
    xyz[2] = 'z';
    err = nc_put_vara_text(NCID, SpatialVID, start, count, xyz);
    if (err != NC_NOERR) {
        CSmallString error;
        error << "unable to set spatial VID 'x', 'y' and 'z' (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    start[0] = 0;
    count[0] = 3;
    xyz[0] = 'a';
    xyz[1] = 'b';
    xyz[2] = 'c';
    err = nc_put_vara_text(NCID, CellSpatialVID, start, count, xyz);
    if (err != NC_NOERR) {
        CSmallString error;
        error << "unable to set spatial cell VID 'a', 'b' and 'c' (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    start[0] = 0;
    start[1] = 0;
    count[0] = 3;
    count[1] = 5;
    err = nc_put_vara_text(NCID, CellAngularVID, start, count, abc);
    if (err != NC_NOERR) {
        CSmallString error;
        error << "unable to set angular cell VID 'alpha', 'beta ' and 'gamma' (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

int CNetCDFTraj::ReadSnapshot(CAmberRestart* p_snap)
{
    if( Mode != AMBER_TRAJ_READ ){
        ES_ERROR("illegal mode, it should be AMBER_TRAJ_READ");
        return(-1);
    }

    if( p_snap == NULL ){
        INVALID_ARGUMENT("p_snap == NULL");
    }

    if( p_snap->GetTopology() == NULL ) {
        ES_ERROR("snapshot does not have assigned topology");
        return(-1);
    }

    if( p_snap->GetNumberOfAtoms() != ActualAtoms ) {
        CSmallString error;
        error << "inconsistent number of atoms, trajectory: " << ActualAtoms;
        error << " topology: " << p_snap->GetNumberOfAtoms();
        ES_ERROR(error);
        return(-1);
    }

    bool has_box = p_snap->GetTopology()->BoxInfo.GetType() != AMBER_BOX_NONE;
    if( has_box != HasBox ) {
        CSmallString error;
        error << "topology and snapshot has different info about box presence";
        ES_ERROR(error);
        return(-1);
    }

    if( CoordinateVID < 0 ) {
        CSmallString error;
        error << "ReadHeader must be called before ReadSnapshot";
        ES_ERROR(error);
        return(-1);
    }

    if( Coordinates == NULL ) {
        CSmallString error;
        error << "Coordinates are NULL";
        ES_ERROR(error);
        return(-1);
    }

    if( CurrentSnapshot >= TotalSnapshots ) return(1); // end of trajectory

    int     err;
    size_t  start[3],count[3];

    // coordinates -------------------------------
    start[0] = CurrentSnapshot;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = ActualAtoms;
    count[2] = 3;

    err = nc_get_vara_float(NCID,CoordinateVID,start,count,Coordinates);
    if( err != NC_NOERR ) {
        CSmallString error;
        error << "unable to get coordinates (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(-1);
    }

    int j=0;
    for(int i=0; i < ActualAtoms; i++) {
        CPoint pos;
        pos.x = Coordinates[j++];
        pos.y = Coordinates[j++];
        pos.z = Coordinates[j++];
        p_snap->SetPosition(i,pos);
    }

    // box ---------------------------------------
    if( HasBox ) {
        start[0] = CurrentSnapshot;
        start[1] = 0;
        start[2] = 0;
        count[0] = 1;
        count[1] = 3;
        count[2] = 0;

        err = nc_get_vara_double(NCID,CellLengthVID,start,count,CellLength);
        if( err != NC_NOERR ) {
            CSmallString error;
            error << "unable to get cell length (" << nc_strerror(err) << ")";
            ES_ERROR(error);
            return(-1);
        }
        CPoint tmp;
        tmp.x =  CellLength[0];
        tmp.y =  CellLength[1];
        tmp.z =  CellLength[2];
        p_snap->SetBox(tmp);

        err = nc_get_vara_double(NCID,CellAngleVID,start,count,CellAngle);
        if( err != NC_NOERR ) {
            CSmallString error;
            error << "unable to get cell length (" << nc_strerror(err) << ")";
            ES_ERROR(error);
            return(-1);
        }

        tmp.x =  CellAngle[0];
        tmp.y =  CellAngle[1];
        tmp.z =  CellAngle[2];
        p_snap->SetAngles(tmp);
    }

    // time --------------------------------------
    if( TimeVID >= 0 ) {
        start[0] = CurrentSnapshot;
        count[0] = 1;
        err = nc_get_vara_float(NCID,TimeVID,start,count,&Time);
        if( err != NC_NOERR ) {
            CSmallString error;
            error << "unable to get time (" << nc_strerror(err) << ")";
            ES_ERROR(error);
            return(-1);
        }
        p_snap->SetTime(Time);
    } else {
        p_snap->SetTime(0.0);
    }

    CurrentSnapshot++;

    return(0);
}

//------------------------------------------------------------------------------

bool CNetCDFTraj::WriteSnapshot(CAmberRestart* p_snap)
{
    if( Mode != AMBER_TRAJ_WRITE ){
        ES_ERROR("illegal mode, it should be AMBER_TRAJ_WRITE");
        return(false);
    }

    if( p_snap == NULL ){
        INVALID_ARGUMENT("p_snap == NULL");
    }

    if( p_snap->GetTopology() == NULL ) {
        ES_ERROR("snapshot does not have assigned topology");
        return(false);
    }

    if( p_snap->GetNumberOfAtoms() != ActualAtoms ) {
        CSmallString error;
        error << "inconsistent number of atoms, trajectory: " << ActualAtoms;
        error << " topology: " << p_snap->GetNumberOfAtoms();
        ES_ERROR(error);
        return(false);
    }

    bool has_box = p_snap->GetTopology()->BoxInfo.GetType() != AMBER_BOX_NONE;
    if( has_box != HasBox ) {
        CSmallString error;
        error << "topology and snapshot has different info about box presence";
        ES_ERROR(error);
        return(false);
    }

    if( CoordinateVID < 0 ) {
        CSmallString error;
        error << "WriteHeader must be called before WriteSnapshot";
        ES_ERROR(error);
        return(false);
    }

    if( Coordinates == NULL ) {
        CSmallString error;
        error << "Coordinates are NULL";
        ES_ERROR(error);
        return(false);
    }

    int j = 0;
    for(int i=0; i < ActualAtoms; i++){
        CPoint pos = p_snap->GetPosition(i);
        Coordinates[j++] = pos.x;
        Coordinates[j++] = pos.y;
        Coordinates[j++] = pos.z;
    }

    int err;

    size_t start[3];
    size_t count[3];

    start[0] = CurrentSnapshot;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = ActualAtoms;
    count[2] = 3;
    err = nc_put_vara_float(NCID, CoordinateVID, start, count,Coordinates);
    if( err != NC_NOERR ){
        CSmallString error;
        error << "unable to write coordinates (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    if( HasBox ) {
        start[0] = CurrentSnapshot;
        start[1] = 0;
        start[2] = 0;
        count[0] = 1;
        count[1] = 3;
        count[2] = 0;

        CPoint box = p_snap->GetBox();
        CellLength[0] = box.x;
        CellLength[1] = box.y;
        CellLength[2] = box.z;

        err = nc_put_vara_double(NCID, CellLengthVID, start, count, CellLength);
        if( err != NC_NOERR ){
            CSmallString error;
            error << "unable to write cell lengths (" << nc_strerror(err) << ")";
            ES_ERROR(error);
            return(false);
        }

        CPoint ang = p_snap->GetAngles();
        CellAngle[0] = ang.x;
        CellAngle[1] = ang.y;
        CellAngle[2] = ang.z;

        err = nc_put_vara_double(NCID, CellAngleVID, start, count, CellAngle);
        if( err != NC_NOERR ){
            CSmallString error;
            error << "unable to write cell angles (" << nc_strerror(err) << ")";
            ES_ERROR(error);
            return(false);
        }
    }

    start[0] = CurrentSnapshot;
    count[0] = 1;
    float time = CurrentSnapshot;
    err = nc_put_vara_float(NCID, TimeVID, start, count, &time);
    if( err != NC_NOERR ){
        CSmallString error;
        error << "unable to write time (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    CurrentSnapshot++;
    TotalSnapshots++;

    nc_sync(NCID);

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



