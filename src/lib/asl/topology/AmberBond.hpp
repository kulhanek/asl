#ifndef AmberBondH
#define AmberBondH
/** \ingroup AmberTopology*/
/*! \file AmberBond.hpp */
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

/// bond

class ASL_PACKAGE CAmberBond {
public:
    /// constructor
    CAmberBond(void);

    /// return index of atom involved in bond
    /*! AMBER abreviation: IBH, IB, IBPER
        index starts from zero and has regular counting order, eg. 0,1,2,... !!
        -1 represents illegal value
    */
    int GetIB(void) const;

    /// set index of atom involved in bond
    void SetIB(int index);

    /// return index of atom involved in bond
    /*! AMBER abreviation: JBH, JB, JBPER
        index starts from zero and has regular counting order, eg. 0,1,2,... !!
        -1 represents illegal value
    */
    int GetJB(void) const;

    /// set index of atom involved in bond
    void SetJB(int index);

    /// return index into parameter arrays BondTypes at lambda=1 (initial system)
    /*! AMBER abreviation: ICBH, ICB, ICBPER - first half
        index starts from zero
        -1 represents illegal value
    */
    int GetICB(void) const;

    /// set index into parameter arrays BondTypes at lambda=1 (initial system)
    void SetICB(int index);

    /// return index into parameter arrays BondTypes at lambda=0 (fully perturbed system)
    /*! AMBER abreviation: ICBPER - second half
        index starts from zero
        -1 represents illegal value
    */
    int GetPCB(void) const;

    /// set index into parameter arrays BondTypes at lambda=0 (fully perturbed system)
    void SetPCB(int index);

// section of private data ----------------------------------------------------
private:
    /// index of atom involved in bond
    /*! AMBER abreviation: IBH, IB, IBPER
        index starts from zero and has regular counting order !!
        -1 represents illegal value
    */
    int IB;
    /// index of atom involved in bond
    /*! AMBER abreviation: JBH, JB, JBPER
        index starts from zero and has regular counting order !!
        -1 represents illegal value
    */
    int JB;
    /// index into parameter arrays BondTypes at lambda=1 (initial system)
    /*! AMBER abreviation: ICBH, ICB, ICBPER - first half
        index starts from zero
        -1 represents illegal value
    */
    int ICB;
    /// index into parameter arrays BondTypes at lambda=0 (fully perturbed system)
    /*! AMBER abreviation: ICBPER - second half
        index starts from zero
        -1 represents illegal value
    */
    int PCB;

    friend class CAmberBondList;
};

//---------------------------------------------------------------------------
#endif
