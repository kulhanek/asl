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
    NCID = -1;
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

    CellLengthVID = -1;
    CellLengthDID = -1;
    CellAngleVID = -1;
    CellAngleDID = -1;

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
    if( NCID >= 0 ) nc_close(NCID);
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
        int err = nc_open(name,NC_NOWRITE,&NCID);
        if( err != NC_NOERR ) {
            CSmallString error;
            error << "unable to open file '" << name << "' for reading (" << nc_strerror(err) << ")";
            ES_ERROR(error);
            return(false);
        }
    }

    if( Mode == AMBER_TRAJ_WRITE ) {
        int err = nc_create(name,NC_64BIT_OFFSET,&NCID);
        if( err != NC_NOERR ) {
            CSmallString error;
            error << "unable to create file '" << name << "' for writing (" << nc_strerror(err) << ")";
            ES_ERROR(error);
            return(false);
        }
    }

    return(true);
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

    // sanity check - time ---------------------------------
    TimeVID = GetVariableID("time");
    if( TimeVID < 0 ) {
        return(false);
    }
    CSmallString unit;
    if( GetVariableAttribute(TimeVID,"units",unit) == false ) {
        return(false);
    }
    if( unit != "picosecond" ) {
        CSmallString error;
        error << "incorrect unit for time (" << unit << "), requested picosecond";
        ES_ERROR(error);
        return(false);
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
    DefineDimension(AMBER_NETCDF_CELL_SPATIAL, 3, &CellLengthDID);
    DefineDimension(AMBER_NETCDF_CELL_ANGULAR, 3, &CellAngleDID);

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

    dimensionID[0] = CellLengthDID;
    DefineVariable(AMBER_NETCDF_CELL_SPATIAL, NC_CHAR, 1, dimensionID, &CellLengthVID);

    dimensionID[0] = CellAngleDID;
    dimensionID[1] = LabelDID;
    DefineVariable(AMBER_NETCDF_CELL_ANGULAR, NC_CHAR, 2, dimensionID, &CellAngleVID);

    // set up box coords
    HasBox = p_top->BoxInfo.GetType() != AMBER_BOX_NONE;

    if( HasBox ) {
        dimensionID[0] = TimeDID;
        dimensionID[1] = CellLengthDID;
        DefineVariable("cell_lengths", NC_DOUBLE, 2, dimensionID, &CellLengthVID);
        PutAttributeText(CellLengthVID, "units", "angstrom");

        dimensionID[1] = CellAngleDID;
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
    err = nc_put_vara_text(NCID, CellLengthVID, start, count, xyz);
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
    err = nc_put_vara_text(NCID, CellAngleVID, start, count, abc);
    if (err != NC_NOERR) {
        CSmallString error;
        error << "unable to set angular cell VID 'alpha', 'beta ' and 'gamma' (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CNetCDFTraj::ReadSnapshot(CAmberRestart* p_snap)
{
    if( Mode != AMBER_TRAJ_READ ){
        ES_ERROR("illegal mode, it should be AMBER_TRAJ_READ");
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
        error << "ReadHeader must be called before ReadSnapshot";
        ES_ERROR(error);
        return(false);
    }

    if( Coordinates == NULL ) {
        CSmallString error;
        error << "Coordinates are NULL";
        ES_ERROR(error);
        return(false);
    }

    if( CurrentSnapshot >= TotalSnapshots ) return(false); // end of trajectory

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
        return(false);
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
            return(false);
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
            return(false);
        }

        tmp.x =  CellAngle[0];
        tmp.y =  CellAngle[1];
        tmp.z =  CellAngle[2];
        p_snap->SetAngles(tmp);
    }

    // time --------------------------------------
    start[0] = CurrentSnapshot;
    count[0] = 1;
    err = nc_get_vara_float(NCID,TimeVID,start,count,&Time);
    if( err != NC_NOERR ) {
        CSmallString error;
        error << "unable to get time (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }
    p_snap->SetTime(Time);

    CurrentSnapshot++;

    return(true);
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

bool CNetCDFTraj::IsNetCDFFile(const CSmallString& name)
{
    int ncid;
    int err = nc_open(name,NC_NOWRITE,&ncid);
    if( err != NC_NOERR ) {
        CSmallString error;
        error << "file '" << name << "' is not NetCDF file (" << nc_strerror(err) << ")";
        ES_TRACE_ERROR(error);
        return(false);
    }
    nc_close(ncid);
    return(true);
}

//------------------------------------------------------------------------------

int CNetCDFTraj::GetDimensionInfo(const char* p_attribute, int* p_length)
{
    int err, dimID;
    size_t slength;

    *p_length = 0;

    err = nc_inq_dimid(NCID, p_attribute, &dimID);
    if (err != NC_NOERR) {
        CSmallString error;
        error << "error on ID of attribute " << p_attribute << " (";
        error << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(-1);
    } else {
        err = nc_inq_dimlen(NCID, dimID, &slength);
        if (err != NC_NOERR) {
            CSmallString error;
            error << "error on value of attribute " << p_attribute << " (";
            error << nc_strerror(err) << ")";
            ES_ERROR(error);
            return(-1);
        }
    }

    *p_length = (int) slength;
    return(dimID);
}

//------------------------------------------------------------------------------

bool CNetCDFTraj::DefineDimension(const char *name, int length, int *dimidp)
{
    int err;
    err = nc_def_dim(NCID, name, length, dimidp);

    if (err != NC_NOERR) {
        CSmallString error;
        error << "error to define dimmension " << name << " (";
        error << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CNetCDFTraj::PutAttributeText(int vid, const char *attribute, const char *text)
{
    int err;
    err = nc_put_att_text(NCID, vid, attribute, strlen(text), text);

    if (err != NC_NOERR) {
        CSmallString error;
        error << "error to put attribute " << attribute << " (";
        error << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }
    return(true);
}

//------------------------------------------------------------------------------

int CNetCDFTraj::GetVariableID(const char* p_variable)
{
    int err;
    int varID = 0;

    err = nc_inq_varid(NCID,p_variable,&varID);
    if (err != NC_NOERR) {
        CSmallString error;
        error << "error on ID of variable " << p_variable << " (";
        error << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(-1);
    }

    return(varID);
}

//------------------------------------------------------------------------------

bool CNetCDFTraj::DefineVariable(const char *name, nc_type xtype, int ndims, int dimids[], int *varidp)
{
    int err;
    err = nc_def_var(NCID, name, xtype, ndims, dimids, varidp);

    if( err != NC_NOERR ) {
        CSmallString error;
        error << "error to define variable " << name << " (";
        error << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }
    return(false);
}

//------------------------------------------------------------------------------

bool CNetCDFTraj::GetVariableAttribute(int varid,const char* p_attribute,CSmallString& text)
{
    int     err;
    int     i;
    size_t  ist;

    text = NULL;

    err = nc_inq_varnatts(NCID,varid,&i);
    if( i <= 0 ) {
        CSmallString error;
        error << "no attributes for variable, varid: " << varid << ", attr: " << p_attribute;
        ES_ERROR(error);
        return(false);
    }

    err = nc_inq_attlen(NCID,varid,p_attribute,&ist);
    if( err != NC_NOERR ) {
        CSmallString error;
        error << "unable to get attribute length of variable (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    text.SetLength(ist);
    nc_get_att_text(NCID,varid,p_attribute,text.GetBuffer());
    if( err != NC_NOERR ) {
        CSmallString error;
        error << "unable to read attribute of variable (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



