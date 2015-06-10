#ifndef AmberResidueH
#define AmberResidueH
/** \ingroup AmberTopology*/
/*! \file AmberResidue.hpp */
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
#include <stdio.h>

//------------------------------------------------------------------------------

class CAmberAtomList;
class CAmberTopology;

//------------------------------------------------------------------------------

/// residue description for topology

class ASL_PACKAGE CAmberResidue {
public:
    CAmberResidue(void);

    /// return residue name
    const char* GetName(bool pert=false) const;

    /// set residue name
    void  SetName(const char* p_name,bool pert=false);

    /// return index of the first atom
    /*! index starts from zero and has regular counting order
    */
    int   GetFirstAtomIndex(void) const;

    /// set index of the first atom
    void  SetFirstAtomIndex(int index);

    /// return number of atoms which belongs to residue
    int   GetNumberOfAtoms(void) const;

    /// return index of residue
    /*! index starts from zero and has regular counting order
    */
    int   GetIndex(void) const;

    /// return number of bonds including hydrogen which belongs to residue
    int   GetNumberOfBondsWithHydrogen(void) const;

    /// return amber topology
    CAmberTopology* GetTopology(void) const;

// section of private data ----------------------------------------------------
private:
    CAmberTopology* Topology;

// AMBER original information ----------------------------------------------
    /// LABRES : the residue labels at lambda=1
    char    LABRES[5];

    /// LABRES : the residue labels at lambda=0
    char    PERRES[5];

    /// IPRES  : atoms in each residue are listed for atom "i" in IPRES(i) to IPRES(i+1)-1
    int     IPRES;

// additional information --------------------------------------------------
    int     Index;
    int     NumOfAtoms;
    int     NumOfBondsWithHydrogen;

    friend class CAmberResidueList;
    friend class CAmberTopology;
};

//------------------------------------------------------------------------------

#endif
