#ifndef AmberMaskResiduesH
#define AmberMaskResiduesH
/** \ingroup AmberMask*/
/*! \file AmberMaskResidues.hpp */
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

#include <SmallString.hpp>
#include <SimpleList.hpp>
#include <ASLMainHeader.hpp>
#include <vector>

//---------------------------------------------------------------------------

class CAmberTopology;
class CAmberRestart;
class CAmberResidue;
class CAmberMaskRSelection;

//---------------------------------------------------------------------------

/// mask support class
/*!
 fully compatible with AMBER 9.0 but only residues can be selected
*/

class ASL_PACKAGE CAmberMaskResidues {
public:
    CAmberMaskResidues(void);
    ~CAmberMaskResidues(void);

    // topology assigment ------------------------------------------------------
    bool            AssignTopology(CAmberTopology* p_top);
    CAmberTopology* GetTopology(void) const;

    // coordinates assigment ---------------------------------------------------
    bool            AssignCoordinates(CAmberRestart* p_crd);
    CAmberRestart*  GetCoordinates(void) const;

    // mask setup --------------------------------------------------------------
    /// select all residues
    bool SelectAllResidues(void);

    /// clear selection
    void Reset(void);

    /// set mask from string specification
    bool                    SetMask(const CSmallString& mask);

    /// set mask from the first line of file maskfile
    bool                    SetMaskFromFile(const CSmallString& maskfile);

    /// set mask from a line read from file p_file
    bool                    SetMaskFromFile(FILE* p_file);

    /// return current mask specification
    const CSmallString&     GetMask(void) const;

    // results -----------------------------------------------------------------
    /// get total number of residues
    int                 GetNumberOfTopologyResidues(void);

    /// get number of selected residues
    int                 GetNumberOfSelectedResidues(void);

    /// get selected residue - index is atom index from 0 to GetNumberOfTopologyResidues(), if residue is not selected NULL is returned
    CAmberResidue*      GetSelectedResidue(int index);

    /// is residue selected - index is atom index from 0 to GetNumberOfTopologyResidues()
    bool                IsResidueSelected(int index);

    /// return selected atom - index is from 0 to GetNumberOfSelectedResidues()
    CAmberResidue*      GetSelectedResidueCondensed(int index);

    /// print mask info
    void                PrintInfo(bool short_info=false,FILE* fout=NULL);

// section of private data ----------------------------------------------------
private:
    CAmberTopology*                 Topology;
    CAmberRestart*                  Coordinates;
    CSmallString                    Mask;
    CAmberMaskRSelection*           Selection;
    std::vector<CAmberResidue*>     SelectedResidues;
};

//---------------------------------------------------------------------------

#endif
