#ifndef AmberBondTypeH
#define AmberBondTypeH
/** \ingroup AmberTopology*/
/*! \file AmberBondType.hpp */
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

/// parametres of bond

class ASL_PACKAGE CAmberBondType {
public:
    CAmberBondType(void);

    /// return force constant
    /*! AMBER abreviation: RK
        unit of force constant: kcal/mol A**2
    */
    double GetRK(void) const;

    /// change force constant
    /*! AMBER abreviation: RK
        unit of force constant: kcal/mol A**2
    */
    void SetRK(double rk);

    /// return equilibrium bond length
    /*! AMBER abreviation: REQ
        unit of length: A
    */
    double GetREQ(void) const;

    /// change equilibrium bond length
    /*! AMBER abreviation: REQ
        unit of length: A
    */
    void SetREQ(double req);

// section of private data ----------------------------------------------------
private:
    /// force constant
    /*! AMBER abreviation: RK
        unit of force constant: kcal/mol A**2
    */
    double  RK;
    /// equilibrium bond length
    /*! AMBER abreviation: REQ
        unit of length: A
    */
    double  REQ;

    friend class CAmberBondList;
};

//------------------------------------------------------------------------------
#endif
