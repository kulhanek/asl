#ifndef AmberMaskAtomsH
#define AmberMaskAtomsH
/** \ingroup AmberMask*/
/*! \file AmberMaskAtoms.hpp */
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
#include <ASLMainHeader.hpp>
#include <vector>

//---------------------------------------------------------------------------

class CAmberTopology;
class CAmberRestart;
class CAmberAtom;
class CAmberMaskASelection;

//---------------------------------------------------------------------------

/// mask support class
/*!
 mask support is now fully compatible with AMBER 9.0
*/

class ASL_PACKAGE CAmberMaskAtoms {
public:
    CAmberMaskAtoms(void);
    ~CAmberMaskAtoms(void);

    // topology assigment ------------------------------------------------------
    bool            AssignTopology(CAmberTopology* p_top);
    CAmberTopology* GetTopology(void) const;

    // coordinates assigment ---------------------------------------------------
    bool            AssignCoordinates(CAmberRestart* p_crd);
    CAmberRestart*  GetCoordinates(void) const;

    // mask setup --------------------------------------------------------------
    /// select all atoms
    bool SelectAllAtoms(void);

    /// clear selection
    void Reset(void);

    /// set mask from string specification
    bool SetMask(const CSmallString& mask);

    /// set mask from the first line of file maskfile
    bool SetMaskFromFile(const CSmallString& maskfile);

    /// set mask from a line read from file p_file
    bool SetMaskFromFile(FILE* p_file);

    /// return current mask specification
    const CSmallString&     GetMask(void) const;

    // results -----------------------------------------------------------------
    /// get total number of atoms
    int             GetNumberOfTopologyAtoms(void);

    /// get number of selected atoms
    int             GetNumberOfSelectedAtoms(void);

    /// get selected atom - index is atom index from 0 to GetNumberOfTopologyAtoms(), if atom is not selected NULL is returned
    CAmberAtom*     GetSelectedAtom(int index);

    /// is atom selected - index is atom index from 0 to GetNumberOfTopologyAtoms()
    bool            IsAtomSelected(int index);

    /// return selected atom - index is from 0 to GetNumberOfSelectedAtoms()
    CAmberAtom*     GetSelectedAtomCondensed(int index);

    /// print mask info
    void            PrintInfo(bool short_info=false,FILE* fout=NULL);

// section of private data ----------------------------------------------------
private:
    CAmberTopology*             Topology;
    CAmberRestart*              Coordinates;
    CSmallString                Mask;
    CAmberMaskASelection*       Selection;
    std::vector<CAmberAtom*>    SelectedAtoms;
};

//---------------------------------------------------------------------------

#endif
