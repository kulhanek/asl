#ifndef AmberAngleH
#define AmberAngleH
/** \ingroup AmberTopology*/
/*! \file AmberAngle.hpp */
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

/// angle description

class ASL_PACKAGE CAmberAngle {
public:
    CAmberAngle(void);

    /// return index of atom involved in angle
    /*! AMBER abreviation: ITH, IT, ITPER
        index starts from zero and has regular counting order, e.g. 0,1,2,... !!
        -1 represents illegal value
    */
    int GetIT(void) const;

    /// set index of atom involved in angle
    void SetIT(int index);

    /// return index of atom involved in angle
    /*! AMBER abreviation: JTH, JT, JTPER
        index starts from zero and has regular counting order, e.g. 0,1,2,... !!
        -1 represents illegal value
    */
    int GetJT(void) const;

    /// set index of atom involved in angle
    void SetJT(int index);

    /// return index of atom involved in angle
    /*! AMBER abreviation: KTH, KT, KTPER
        index starts from zero and has regular counting order, e.g. 0,1,2,... !!
        -1 represents illegal value
    */
    int GetKT(void) const;

    /// set index of atom involved in angle
    void SetKT(int index);

    /// return index into parameter arrays BondTypes at lambda=1 (initial system)
    /*! AMBER abreviation: ICTH, ICT, ICBPER - first half
        index starts from zero
        -1 represents illegal value
    */
    int GetICT(void) const;

    /// set index into parameter arrays BondTypes at lambda=1 (initial system)
    void SetICT(int index);

    /// return index into parameter arrays BondTypes at lambda=0 (fully perturbed system)
    /*! AMBER abreviation: ICBPER - second half
        index starts from zero
        -1 represents illegal value
    */
    int GetPCT(void) const;

    /// set index into parameter arrays BondTypes at lambda=0 (fully perturbed system)
    void SetPCT(int index);

// section of private data ----------------------------------------------------
private:
    /// index of atom involved in angle
    /*! AMBER abreviation: ITH, IT, ITPER
        index starts from zero and has regular counting order !!
        -1 represents illegal value
    */
    int IT;

    /// index of atom involved in angle
    /*! AMBER abreviation: JTH, JT, JTPER
        index starts from zero and has regular counting order !!
        -1 represents illegal value
    */
    int JT;

    /// index of atom involved in angle
    /*! AMBER abreviation: KTH, KT, KTPER
        index starts from zero and has regular counting order !!
        -1 represents illegal value
    */
    int KT;

    /// index into parameter arrays BondTypes at lambda=1 (initial system)
    /*! AMBER abreviation: ICTH, ICT, ICBPER - first half
        index starts from zero
        -1 represents illegal value
    */
    int ICT;

    /// index into parameter arrays BondTypes at lambda=0 (fully perturbed system)
    /*! AMBER abreviation: ICBPER - second half
        index starts from zero
        -1 represents illegal value
    */
    int PCT;

private:
    CAmberAngle(const CAmberAngle& copy);

    friend class CAmberAngleList;
};

//---------------------------------------------------------------------------

#endif
