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

#include <NetCDFFile.hpp>
#include <ErrorSystem.hpp>
#include <string.h>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CNetCDFFile::CNetCDFFile(void)
{
    NCID = -1;
}

//---------------------------------------------------------------------------

CNetCDFFile::~CNetCDFFile(void)
{
    if( NCID >= 0 ) nc_close(NCID);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CNetCDFFile::Open(const CSmallString& name,char mode)
{
    if( NCID >= 0 ) {
        ES_ERROR("file is already opened");
        return(false);
    }

    if( mode == 'r' ) {
        int err = nc_open(name,NC_NOWRITE,&NCID);
        if( err != NC_NOERR ) {
            CSmallString error;
            error << "unable to open file '" << name << "' for reading (" << nc_strerror(err) << ")";
            ES_ERROR(error);
            return(false);
        }
        return(true);
    }

    if( mode == 'w' ) {
        int err = nc_create(name,NC_64BIT_OFFSET,&NCID);
        if( err != NC_NOERR ) {
            CSmallString error;
            error << "unable to create file '" << name << "' for writing (" << nc_strerror(err) << ")";
            ES_ERROR(error);
            return(false);
        }
        return(true);
    }

    ES_ERROR("unsupported mode");
    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CNetCDFFile::IsNetCDFFile(const CSmallString& name)
{
    int ncid;
    int err = nc_open(name,NC_NOWRITE,&ncid);
    if( err != NC_NOERR ) {
//        CSmallString error;
//        error << "file '" << name << "' is not NetCDF file (" << nc_strerror(err) << ")";
//        ES_TRACE_ERROR(error);
        return(false);
    }
    nc_close(ncid);
    return(true);
}

//------------------------------------------------------------------------------

int CNetCDFFile::GetDimensionInfo(const char* p_attribute, int* p_length)
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

bool CNetCDFFile::DefineDimension(const char *name, int length, int *dimidp)
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

bool CNetCDFFile::PutAttributeText(int vid, const char *attribute, const char *text)
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

bool CNetCDFFile::PutAttributeValue(int vid, const char *attribute, double value)
{
    int err;
    err = nc_put_att_double(NCID, vid, attribute,NC_DOUBLE,1, &value);

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

int CNetCDFFile::GetVariableID(const char* p_variable)
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

bool CNetCDFFile::DefineVariable(const char *name, nc_type xtype, int ndims, int dimids[], int *varidp)
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

bool CNetCDFFile::GetVariableAttribute(int varid,const char* p_attribute,CSmallString& text)
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



