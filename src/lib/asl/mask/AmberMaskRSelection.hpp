#ifndef AmberMaskRSelectionH
#define AmberMaskRSelectionH
/** \ingroup AmberMask*/
/*! \file AmberMaskRSelection.hpp */
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
#include "maskparser/AmberMaskParser.hpp"

//---------------------------------------------------------------------------

class CAmberMaskResidues;
class CAmberResidue;

//---------------------------------------------------------------------------

/// selection of atoms

class ASL_PACKAGE CAmberMaskRSelection {
public:
    CAmberMaskRSelection(CAmberMaskResidues* p_owner);
    ~CAmberMaskRSelection(void);

    // executive methods  ------------------------------------------------------
    bool ExpandAndReduceTree(struct SExpression* p_expr);

    // results -----------------------------------------------------------------
    int                 GetNumberOfSelectedResidues(void);
    CAmberResidue*      GetSelectedResidue(int index);

// section of private data ----------------------------------------------------
private:
    CAmberMaskResidues*        Owner;
    CAmberResidue**            Residues;

    static bool ExpandAndReduceTree(CAmberMaskRSelection* p_root,
            struct SExpression* p_expr);

    // individual selections
    bool Select(struct SSelection* p_sel);
    void SelectResidueByIndex(int index,int length);
    void SelectResidueByName(const char* p_name);

    bool SelectResidueByDistanceFromOrigin(SOperator dist_oper,double dist);
    bool SelectResidueByDistanceFromCentreOfBox(SOperator dist_oper,double dist);
    bool SelectResidueByDistanceFromList(CAmberMaskRSelection* p_left,
            SOperator dist_oper,double dist);
    bool SelectResidueByDistanceFromCOM(CAmberMaskRSelection* p_left,
            SOperator dist_oper,double dist);
    bool SelectResidueByDistanceFromPlane(CAmberMaskRSelection* p_left,
            SOperator dist_oper,double dist);

};

//---------------------------------------------------------------------------

#endif
