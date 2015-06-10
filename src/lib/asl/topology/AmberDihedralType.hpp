#ifndef AmberDihedralTypeH
#define AmberDihedralTypeH
/** \ingroup AmberTopology*/
/*! \file AmberDihedralType.hpp */
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

/// dihedral type description

class ASL_PACKAGE CAmberDihedralType {
public:
    CAmberDihedralType(void);

    /// return force constant
    /*! AMBER abreviation: PK
        unit of force constant: kcal/mol
    */
    double GetPK(void) const;

    /// change force constant
    /*! AMBER abreviation: PK
        unit of force constant: kcal/mol
    */
    void SetPK(double pk);

    /// return periodicity of dihedral
    /*! AMBER abreviation: PN
        unit: 1
    */
    double GetPN(void) const;

    /// change periodicity of dihedral
    /*! AMBER abreviation: PN
        unit: 1
    */
    void SetPN(double pn);

    /// return phase of dihedral
    /*! AMBER abreviation: PHASE
        unit: rad
    */
    double GetPHASE(void) const;

    /// change phase of dihedral
    /*! AMBER abreviation: PHASE
        unit: rad
    */
    void SetPHASE(double phase);

    /// return 1,4-electrostatic scaling factor
    /*! AMBER abreviation: SCEE_SCALE_FACTOR
        unit: 1
    */
    double GetSCEE(void) const;

    /// change 1,4-electrostatic scaling factor
    /*! AMBER abreviation: SCEE_SCALE_FACTOR
        unit: 1
    */
    void SetSCEE(double scee);

    /// return 1,4-vdW scaling factor
    /*! AMBER abreviation: SCEE_SCNB_FACTOR
        unit: 1
    */
    double GetSCNB(void) const;

    /// change 1,4-vdW scaling factor
    /*! AMBER abreviation: SCEE_SCNB_FACTOR
        unit: 1
    */
    void SetSCNB(double scnb);

// section of private data ----------------------------------------------------
private:
    /// force constant
    /*! AMBER abreviation: PK
        unit of force constant: kcal/mol
    */
    double  PK;

    /// periodicity of dihedral
    /*! AMBER abreviation: PN
        unit: 1
    */
    double  PN;

    /// phase of dihedral
    /*! AMBER abreviation: PHASE
        unit: rad
    */
    double PHASE;

    /// 1,4-electrostatic scaling factor
    /*! AMBER abreviation: SCEE_SCALE_FACTOR
        unit: 1
    */
    double SCEE_SCALE;

    /// 1,4-vdW scaling factor
    /*! AMBER abreviation: SCNB_SCALE_FACTOR
        unit: 1
    */
    double SCNB_SCALE;

    friend class CAmberDihedralList;
};

//---------------------------------------------------------------------------

#endif
