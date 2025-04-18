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
//     This program is distributed","the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with this program; if not, write to the Free Software Foundation, Inc.,
//     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
// =============================================================================

#include <string.h>
#include <stdlib.h>
#include <AmberCMAPList.hpp>
#include <FortranIO.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberCMAPList::CAmberCMAPList(void)
{
    cmap_loaded = false;
    cmap_term_count = 0;
    cmap_type_count = 0;
    cmaps = NULL;
}

//------------------------------------------------------------------------------

CAmberCMAPList::~CAmberCMAPList(void)
{
    FreeFields();
}

//------------------------------------------------------------------------------

void CAmberCMAPList::FreeFields(void)
{
    cmap_loaded = false;
    cmap_term_count = 0;
    cmap_type_count = 0;
    if( cmaps != NULL ) delete[] cmaps;
    cmaps = NULL;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberCMAPList::IsCMAPSection(const char* p_section)
{
    if( strcmp(p_section,"%FLAG CMAP_COUNT") == 0 ) return(true);
    if( strcmp(p_section,"%FLAG CMAP_RESOLUTION") == 0 ) return(true);
    if( strstr(p_section,"%FLAG CMAP_PARAMETER_") != NULL ) return(true);
    if( strcmp(p_section,"%FLAG CMAP_INDEX") == 0 ) return(true);
    return(false);
}

//------------------------------------------------------------------------------

bool CAmberCMAPList::LoadCMAPSection(FILE* p_file,const char* p_section)
{
    CFortranIO fortranio(p_file);

//-----------------------------------
    if( strcmp(p_section,"%FLAG CMAP_COUNT") == 0 ){
        fCMAP_COUNT = fortranio.LoadFormatOfSection();
        fortranio.SetFormat(fCMAP_COUNT);

        if( fortranio.ReadInt(cmap_term_count) == false ) {
            ES_ERROR("unable to load cmap_term_count item");
            return(false);
        }
        if( fortranio.ReadInt(cmap_type_count) == false ) {
            ES_ERROR("unable to load cmap_type_count item");
            return(false);
        }

        if( cmaps != NULL ) delete[] cmaps;
        cmaps = new CAmberCMAP[cmap_type_count];

        cmap_loaded = true;
        return(true);
    }

//-----------------------------------
    if( strcmp(p_section,"%FLAG CMAP_RESOLUTION") == 0 ){
        fCMAP_RESOLUTION = fortranio.LoadFormatOfSection();
        fortranio.SetFormat(fCMAP_RESOLUTION);

        int resolution_length = cmap_type_count;
        if( resolution_length <= 0 ){
            ES_ERROR("illegal resolution_length");
            return(false);
        }

        cmap_resolution.CreateVector(resolution_length);
        for(int i=0; i<resolution_length; i++) {
            if( fortranio.ReadInt(cmap_resolution[i]) == false ) {
                ES_ERROR("unable to load cmap_resolution item");
                return(false);
            }
        }
        cmap_loaded = true;
        return(true);
    }

//-----------------------------------
    if( strstr(p_section,"%FLAG CMAP_PARAMETER_") != NULL ){

        if( cmaps == NULL ){
            ES_ERROR("no cmaps array initialized");
            return(false);
        }

        // find cmap index
        int index = 0;
        for(int i=1; i <= cmap_type_count; i++){
            CSmallString sname = "%FLAG CMAP_PARAMETER_";
            CSmallString sidx;
            sidx.IntToStr(i,"%02d");
            sname << sidx;

            if( strcmp(p_section,sname) == 0) {
                index = i;
                break;
            }
        }

        if( index == 0 ){
            ES_ERROR("cmap index out-of-range");
            return(false);
        }
        index--;

        CAmberCMAP* p_cmap = &cmaps[index];

        CSmallString comment,cmap;
        comment.ReadStringFromFile(p_file);
        p_cmap->cmap_title.ReadStringFromFile(p_file);
        cmap.ReadStringFromFile(p_file);

        p_cmap->fCMAP_PARAMETER = fortranio.LoadFormatOfSection();

        fortranio.SetFormat(p_cmap->fCMAP_PARAMETER);

        int cmap_length = cmap_resolution[index] * cmap_resolution[index];
        if( cmap_length <= 0 ){
            ES_ERROR("illegal cmap_length");
            return(false);
        }

        p_cmap->cmap_data.CreateVector(cmap_length);
        for(int i=0; i<cmap_length; i++) {
            if( fortranio.ReadReal(p_cmap->cmap_data[i]) == false ) {
                ES_ERROR("unable to load cmap_data item");
                return(false);
            }
        }

        cmap_loaded = true;
        return(true);
    }

//-----------------------------------
    if( strcmp(p_section,"%FLAG CMAP_INDEX") == 0 ){
        fCMAP_INDEX = fortranio.LoadFormatOfSection();
        fortranio.SetFormat(fCMAP_INDEX);

        int index_length = cmap_term_count*6;
        if( index_length <= 0 ){
            ES_ERROR("illegal index_length");
            return(false);
        }

        cmap_index.CreateVector(index_length);
        for(int i=0; i<index_length; i++) {
            if( fortranio.ReadInt(cmap_index[i]) == false ) {
                ES_ERROR("unable to load cmap_index item");
                return(false);
            }
        }
        cmap_loaded = true;
        return(true);
    }

   return(false);
}

//------------------------------------------------------------------------------

bool CAmberCMAPList::SaveCMAPSections(FILE* p_file)
{
    if( cmap_loaded == false ) return(true);

    CFortranIO fortranio(p_file);

//-----------------------------------
    if( SaveSectionHeader(p_file,"CMAP_COUNT",fCMAP_COUNT) == false ) return(false);
    fortranio.SetFormat(fCMAP_COUNT);

    if( fortranio.WriteInt(cmap_term_count) == false ) {
        ES_ERROR("unable to save cmap_term_count item");
        return(false);
    }
    if( fortranio.WriteInt(cmap_type_count) == false ) {
        ES_ERROR("unable to save cmap_type_count item");
        return(false);
    }
    fortranio.WriteEndOfSection();

//-----------------------------------
    if( SaveSectionHeader(p_file,"CMAP_RESOLUTION",fCMAP_RESOLUTION) == false ) return(false);

    fortranio.SetFormat(fCMAP_RESOLUTION);
    for(int i=0; i < cmap_type_count; i++){
        if( fortranio.WriteInt(cmap_resolution[i]) == false ) {
            ES_ERROR("unable to save cmap_resolution item");
            return(false);
        }
    }
    fortranio.WriteEndOfSection();

//-----------------------------------
    for(int i=0; i < cmap_type_count; i++){
        CAmberCMAP* p_cmap = &cmaps[i];
        CSmallString sname = "CMAP_PARAMETER_";
        CSmallString sidx;
        sidx.IntToStr(i+1,"%02d");

        sname << sidx;

        if( SaveSectionHeader(p_file,sname,p_cmap->fCMAP_PARAMETER,p_cmap->cmap_title) == false ) return(false);

        int cmap_length = cmap_resolution[i] * cmap_resolution[i];

        fortranio.SetFormat(p_cmap->fCMAP_PARAMETER);
        for(int j=0; j < cmap_length; j++){
            if( fortranio.WriteReal(p_cmap->cmap_data[j]) == false ) {
                ES_ERROR("unable to save cmap_data item");
                return(false);
            }
        }
        fortranio.WriteEndOfSection();
    }

//-----------------------------------
    if( SaveSectionHeader(p_file,"CMAP_INDEX",fCMAP_INDEX) == false ) return(false);

    int index_length = cmap_term_count * 6;

    fortranio.SetFormat(fCMAP_INDEX);
    for(int i=0; i < index_length; i++){
        if( fortranio.WriteInt(cmap_index[i]) == false ) {
            ES_ERROR("unable to save cmap_index item");
            return(false);
        }
    }
    fortranio.WriteEndOfSection();

//-----------------------------------

   return(true);
}

//------------------------------------------------------------------------------

bool CAmberCMAPList::SaveSectionHeader(FILE* p_top,const char* p_section_name,
        const char* p_section_format)
{
    int outputlen;

    if( (outputlen = fprintf(p_top,"%%FLAG %s",p_section_name)) <= 0 ) {
        CSmallString    error;
        error << "unable write header of %%FLAG " << p_section_name << " section";
        ES_ERROR(error);
        return(false);
    }

    // for(int i = outputlen; i < 80; i++) fputc(' ',p_top);
    fputc('\n',p_top);

    if( (outputlen = fprintf(p_top,"%%FORMAT(%s)",(char*)p_section_format)) <= 0 ) {
        CSmallString    error;
        error << "unable write format of %%FLAG " << p_section_name << " section";
        ES_ERROR(error);
        return(false);
    }

    // for(int i = outputlen; i < 80; i++) fputc(' ',p_top);
    fputc('\n',p_top);

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberCMAPList::SaveSectionHeader(FILE* p_top,const char* p_section_name,
        const char* p_section_format,const char* p_comment)
{
    int outputlen;

    if( (outputlen = fprintf(p_top,"%%FLAG %s",p_section_name)) <= 0 ) {
        CSmallString    error;
        error << "unable write header of %%FLAG " << p_section_name << " section";
        ES_ERROR(error);
        return(false);
    }

    // for(int i = outputlen; i < 80; i++) fputc(' ',p_top);
    fputc('\n',p_top);

    if( (outputlen = fprintf(p_top,"%%COMMENT %s CMAP",p_comment)) <= 0 ) {
        CSmallString    error;
        error << "unable write header of %%COMMENT " << p_comment << " CMAP section";
        ES_ERROR(error);
        return(false);
    }

    // for(int i = outputlen; i < 80; i++) fputc(' ',p_top);
    fputc('\n',p_top);

    if( (outputlen = fprintf(p_top,"%%FORMAT(%s)",(char*)p_section_format)) <= 0 ) {
        CSmallString    error;
        error << "unable write format of %%FLAG " << p_section_name << " section";
        ES_ERROR(error);
        return(false);
    }

    // for(int i = outputlen; i < 80; i++) fputc(' ',p_top);
    fputc('\n',p_top);

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberCMAPList::operator = (const CAmberCMAPList& src)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




