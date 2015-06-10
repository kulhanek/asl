#ifndef AmberAngleTypeH
#define AmberAngleTypeH
/** \ingroup AmberTopology*/
/*! \file AmberAngleType.hpp */
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

//---------------------------------------------------------------------------
//###########################################################################
//---------------------------------------------------------------------------

/// angle type

class ASL_PACKAGE CAmberAngleType {
public:
    CAmberAngleType(void);

    /// return force constant
    /*! AMBER abreviation: TK
        unit of force constant: kcal/mol
    */
    double GetTK(void) const;

    /// change force constant
    /*! AMBER abreviation: TK
        unit of force constant: kcal/mol
    */
    void SetTK(double tk);

    /// return equilibrium angle
    /*! AMBER abreviation: TEQ
        unit of equilibrium angle: rad
    */
    double GetTEQ(void) const;

    /// change equilibrium angle
    /*! AMBER abreviation: TEQ
        unit of equilibrium angle: rad
    */
    void SetTEQ(double teq);

// section of private data ----------------------------------------------------
private:
    /// force constant
    /*! AMBER abreviation: TK
        unit of force constant: kcal/mol
    */
    double  TK;

    /// equilibrium angle
    /*! AMBER abreviation: TEQ
        unit of equilibrium angle: rad
    */
    double  TEQ;

    friend class CAmberAngleList;
};

//---------------------------------------------------------------------------

#endif
