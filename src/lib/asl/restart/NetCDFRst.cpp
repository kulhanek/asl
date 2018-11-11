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

#include <NetCDFRst.hpp>
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

CNetCDFRst::CNetCDFRst(void)
{
    Mode = 'r';

    SpatialVID = -1;
    SpatialDID = -1;

    Spatial = 0;
    ActualAtoms = 0;

    CoordinateVID = -1;
    CoordinateDID = -1;
    Coordinates = NULL;

    VelocityVID = -1;
    VelocityDID = -1;
    Velocities = NULL;

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

CNetCDFRst::~CNetCDFRst(void)
{
    if( Coordinates != NULL ) delete[] Coordinates;
    if( Velocities != NULL ) delete[] Velocities;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CNetCDFRst::Open(const CSmallString& name,char mode)
{
    if( NCID >= 0 ) {
        ES_ERROR("file is already opened");
        return(false);
    }
    Mode = mode;
    return(CNetCDFFile::Open(name,mode));
}

//------------------------------------------------------------------------------

bool CNetCDFRst::ReadHeader(CAmberTopology* p_top)
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

    if( Conventions != "AMBERRESTART" ) {
        CSmallString error;
        error << "illegal conventions '" << Conventions << "', expecting 'AMBERRESTART'";
        ES_ERROR(error);
        return(false);
    }
    if( ConventionVersion != "1.0" ) {
        CSmallString error;
        error << "illegal convention version '" << ConventionVersion << "', expecting '1.0'";
        ES_ERROR(error);
        return(false);
    }

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

    if( Coordinates != NULL ) {
        delete Coordinates;
    }
    Coordinates = new double[ActualAtoms*3];

    if( Velocities != NULL ) {
        delete Velocities;
    }
    Velocities = new double[ActualAtoms*3];

    // sanity check - box ----------------------------------
    if( HasBox ){
        // optional
        CellLengthVID = GetVariableID("cell_lengths");
        CellAngleVID = GetVariableID("cell_angles");
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CNetCDFRst::WriteHeader(CAmberTopology* p_top)
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

    if( Coordinates != NULL ) {
        delete Coordinates;
    }
    Coordinates = new double[ActualAtoms*3];

    if( Velocities != NULL ) {
        delete Velocities;
    }
    Velocities = new double[ActualAtoms*3];

    // global dimmensions
    DefineDimension(AMBER_NETCDF_FRAME, NC_UNLIMITED, &TimeDID);
    DefineDimension(AMBER_NETCDF_SPATIAL, 3, &SpatialDID);
    DefineDimension(AMBER_NETCDF_ATOM, ActualAtoms, &CoordinateDID);
    DefineDimension(AMBER_NETCDF_ATOM, ActualAtoms, &VelocityDID);
    DefineDimension(AMBER_NETCDF_LABEL, AMBER_NETCDF_LABELLEN, &LabelDID);
    DefineDimension(AMBER_NETCDF_CELL_SPATIAL, 3, &CellSpatialDID);
    DefineDimension(AMBER_NETCDF_CELL_ANGULAR, 3, &CellAngularDID);

    // put global attributes
    Conventions =  "AMBERRESTART";
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

    dimensionID[0] = CoordinateDID;
    dimensionID[1] = SpatialDID;
    DefineVariable(AMBER_NETCDF_COORDS, NC_FLOAT, 2, dimensionID, &CoordinateVID);
    PutAttributeText(CoordinateVID, "units", "angstrom");

    dimensionID[0] = VelocityDID;
    dimensionID[1] = SpatialDID;
    DefineVariable(AMBER_NETCDF_COORDS, NC_FLOAT, 2, dimensionID, &VelocityVID);
    PutAttributeText(VelocityVID, "units", "angstrom/picosecond");
// velocities:scale_factor = 20.455 ; FIXME

    dimensionID[0] = CellSpatialDID;
    DefineVariable(AMBER_NETCDF_CELL_SPATIAL, NC_CHAR, 1, dimensionID, &CellSpatialVID);

    dimensionID[0] = CellAngularDID;
    dimensionID[1] = LabelDID;
    DefineVariable(AMBER_NETCDF_CELL_ANGULAR, NC_CHAR, 2, dimensionID, &CellAngularVID);

    // set up box coords
    HasBox = p_top->BoxInfo.GetType() != AMBER_BOX_NONE;

    if( HasBox ) {
        dimensionID[0] = CellSpatialDID;
        DefineVariable("cell_lengths", NC_DOUBLE, 1, dimensionID, &CellLengthVID);
        PutAttributeText(CellLengthVID, "units", "angstrom");

        dimensionID[0] = CellAngularDID;
        DefineVariable("cell_angles", NC_DOUBLE, 1, dimensionID, &CellAngleVID);
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

bool CNetCDFRst::ReadSnapshot(CAmberRestart* p_snap)
{
    if( Mode != 'r' ){
        ES_ERROR("illegal mode, it should be 'r'");
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

    int     err;
    size_t  start[2],count[2];

    // coordinates -------------------------------
    start[0] = 0;
    count[0] = ActualAtoms;
    start[1] = 0;
    count[1] = 3;

    err = nc_get_vara_double(NCID,CoordinateVID,start,count,Coordinates);
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

    // velocities -------------------------------
    err = nc_get_vara_double(NCID,VelocityVID,start,count,Velocities);
    if( err != NC_NOERR ) {
        CSmallString error;
        error << "unable to get velocities (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    j=0;
    for(int i=0; i < ActualAtoms; i++) {
        CPoint vel;
        vel.x = Velocities[j++];
        vel.y = Velocities[j++];
        vel.z = Velocities[j++];
        p_snap->SetVelocity(i,vel);
    }

    // box ---------------------------------------
    if( HasBox ) {
        start[0] = 0;
        count[0] = 3;

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
    start[0] = 0;
    count[0] = 1;
    err = nc_get_vara_double(NCID,TimeVID,start,count,&Time);
    if( err != NC_NOERR ) {
        CSmallString error;
        error << "unable to get time (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }
    p_snap->SetTime(Time);

    return(true);
}

//------------------------------------------------------------------------------

bool CNetCDFRst::WriteSnapshot(CAmberRestart* p_snap)
{
    if( Mode != 'w' ){
        ES_ERROR("illegal mode, it should be 'w'");
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

    if( Velocities == NULL ) {
        CSmallString error;
        error << "Velocities are NULL";
        ES_ERROR(error);
        return(false);
    }

    j = 0;
    for(int i=0; i < ActualAtoms; i++){
        CPoint vel = p_snap->GetVelocity(i);
        Velocities[j++] = vel.x;
        Velocities[j++] = vel.y;
        Velocities[j++] = vel.z;
    }

    int err;

    size_t start[2];
    size_t count[2];

    start[0] = 0;
    count[0] = ActualAtoms;
    start[1] = 0;
    count[1] = 3;
    err = nc_put_vara_double(NCID, CoordinateVID, start, count,Coordinates);
    if( err != NC_NOERR ){
        CSmallString error;
        error << "unable to write coordinates (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    err = nc_put_vara_double(NCID, VelocityVID, start, count,Velocities);
    if( err != NC_NOERR ){
        CSmallString error;
        error << "unable to write velocities (" << nc_strerror(err) << ")";
        ES_ERROR(error);
        return(false);
    }

    if( HasBox ) {
        start[0] = 0;
        count[0] = 3;
        start[1] = 0;
        count[1] = 0;

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

    nc_sync(NCID);

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



