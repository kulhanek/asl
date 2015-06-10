#ifndef AmberResidueListH
#define AmberResidueListH
/** \ingroup AmberTopology*/
/*! \file AmberResidueList.hpp */
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
#include <AmberResidue.hpp>
#include <AmberAtomList.hpp>

//------------------------------------------------------------------------------

class CAmberTopology;

//------------------------------------------------------------------------------

/// residue list in topology

class ASL_PACKAGE CAmberResidueList {
public:
    CAmberResidueList(CAmberTopology* p_top);
    ~CAmberResidueList(void);

    /// init all fields - old data are destroyed
    void InitFields(int iNRES, int iNMXRS);

    /// prepare for new data
    void FreeFields(void);

    /// returns number of residues
    int            GetNumberOfResidues(void) const;

    /// return pointer to residue
    /*! it cannot return NULL if index is in legal range
    */
    CAmberResidue* GetResidue(int index) const;

// sorted residue list
    /// prepare sorted residue list
    bool PrepareSortedResidues(void);

    /// free memory associated with sorted list
    bool FreeSortedResidueList(void);

    /// return pointer to residue from sorted list
    /*! return NULL if sorted list is not prepared
    */
    CAmberResidue* GetSortedResidue(int index) const;

    /// reinit atom residue pointers
    void ReinitAtomResiduePointers(CAmberAtomList* p_atomlist);

    /// overload assigment operator
    void operator = (const CAmberResidueList& src);

// section of private data ----------------------------------------------------
private:
    CAmberTopology* Topology;

    /// NRES  :number of residues
    int             NRES;

    /// NMXRS  : number of atoms in the largest residue
    int             NMXRS;

    CAmberResidue*  Residues;
    CAmberResidue** SortedResidues;

    bool LoadResidueNames(FILE* p_file,const char* p_format);
    bool LoadResidueIPRES(FILE* p_file,CAmberAtomList* p_atomlist,const char* p_format);
    bool LoadResiduePertNames(FILE* p_file,const char* p_format);

    bool SaveResidueNames(FILE* p_file,const char* p_format);
    bool SaveResidueIPRES(FILE* p_file,CAmberAtomList* p_atomlist,const char* p_format);
    bool SaveResiduePertNames(FILE* p_file,const char* p_format);

    static int CompareItem(const void *p_left, const void *p_right);

    friend class CAmberTopology;
};

//---------------------------------------------------------------------------

#endif
