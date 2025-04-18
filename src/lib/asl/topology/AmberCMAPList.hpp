#ifndef AmberCMAPListH
#define AmberCMAPListH
/** \ingroup AmberTopology*/
/*! \file AmberCMAPList.hpp */
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

#include <stdio.h>
#include <ASLMainHeader.hpp>
#include <SmallString.hpp>
#include <SimpleVector.hpp>
#include <AmberCMAPList.hpp>

//---------------------------------------------------------------------------

class ASL_PACKAGE CAmberCMAP {
public:
    CSmallString            fCMAP_PARAMETER;
    CSmallString            cmap_title;
    CSimpleVector<double>   cmap_data;
};

//---------------------------------------------------------------------------

/// CMAP list for topology

class ASL_PACKAGE CAmberCMAPList {
public:
    CAmberCMAPList(void);
    ~CAmberCMAPList(void);

    /// prepare for new data
    void FreeFields(void);

    /// overload assigment operator
    void operator = (const CAmberCMAPList& src);

// section of private data ----------------------------------------------------
private:
    bool    cmap_loaded;

    // CMAP_COUNT
    int     cmap_term_count;
    int     cmap_type_count;

    // CMAP_RESOLUTION
    CSimpleVector<int>  cmap_resolution;

    // CMAP_INDEX
    CSimpleVector<int>  cmap_index;

    // CMAPS
    CAmberCMAP*         cmaps;

    // formats
    CSmallString fCMAP_COUNT;
    CSmallString fCMAP_RESOLUTION;
    CSmallString fCMAP_INDEX;

    bool IsCMAPSection(const char* p_section);
    bool LoadCMAPSection(FILE* p_file,const char* p_section);
    bool SaveCMAPSections(FILE* p_file);

    bool SaveSectionHeader(FILE* p_top,const char* p_section_name,
            const char* p_section_format);
    bool SaveSectionHeader(FILE* p_top,const char* p_section_name,
            const char* p_section_format,const char* p_comment);

    friend class CAmberTopology;
    friend class CAmberSubTopology;
};

//---------------------------------------------------------------------------

#endif
