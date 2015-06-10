#ifndef AmberDihedralH
#define AmberDihedralH
/** \ingroup AmberTopology*/
/*! \file AmberDihedral.hpp */
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

/// dihedral description

class ASL_PACKAGE CAmberDihedral {
public:
    CAmberDihedral(void);

    /// return index of atom involved in dihedral angle
    /*! AMBER abreviation: IPH, IP, IPPER
        index starts from zero and has regular counting order, e.g. 0,1,2,... !!
        -1 represents illegal value
    */
    int GetIP(void) const;

    /// set index of atom involved in dihedral angle
    void SetIP(int index);

    /// return index of atom involved in dihedral angle
    /*! AMBER abreviation: JPH, JP, JPPER
        index starts from zero and has regular counting order, e.g. 0,1,2,... !!
        -1 represents illegal value
    */
    int GetJP(void) const;

    /// set index of atom involved in dihedral angle
    void SetJP(int index);

    /// return index of atom involved in dihedral angle
    /*! AMBER abreviation: KPH, KP, KPPER
        AMBER meaning for KP < 0 is substituted with Type = -1
        index starts from zero and has regular counting order, e.g. 0,1,2,... !!
        -1 represents illegal value
    */
    int GetKP(void) const;

    /// set index of atom involved in dihedral angle
    void SetKP(int index);

    /// return index of atom involved in dihedral angle
    /*! AMBER abreviation: LPH, LP, LPPER
        AMBER meaning for LP < 0 is substituted with Type = 1
        index starts from zero and has regular counting order, e.g. 0,1,2,... !!
        -1 represents illegal value
    */
    int GetLP(void) const;

    /// return index of atom involved in dihedral angle
    void SetLP(int index);

    /// return index into parameter arrays DihedralTypes at lambda=1 (initial system)
    /*! AMBER abreviation: ICPH, ICP, ICPPER - first half
        index starts from zero
        -1 represents illegal value
    */
    int GetICP(void) const;

    /// set index into parameter arrays DihedralTypes at lambda=1 (initial system)
    void SetICP(int index);

    /// return index into parameter arrays DihedralTypes at lambda=0 (fully perturbed system)
    /*! AMBER abreviation: ICPPER - second half
        index starts from zero
        -1 represents illegal value
    */
    int GetPCP(void) const;

    /// set index into parameter arrays DihedralTypes at lambda=0 (fully perturbed system)
    void SetPCP(int index);

    /// return type of dihedral angle
    /*! 0 - normal dihedral angle
        1 - improper dihedral angle
       -1 - terminal dihedral angle
       -2 - terminal and improper dihedral angle
    */
    int GetType(void) const;

    /// set type of dihedral angle
    void SetType(int type);

// section of private data ----------------------------------------------------
private:
    /// index of atom involved in dihedral angle
    /*! AMBER abreviation: IPH, IP, IPPER
        index starts from zero and has regular counting order !!
        -1 represents illegal value
    */
    int     IP;

    /// index of atom involved in dihedral angle
    /*! AMBER abreviation: JPH, JP, JPPER
        index starts from zero and has regular counting order !!
        -1 represents illegal value
    */
    int     JP;

    /// index of atom involved in dihedral angle
    /*! AMBER abreviation: KPH, KP, KPPER
        AMBER meaning for KP < 0 is substituted with Type = -1
        index starts from zero and has regular counting order !!
        -1 represents illegal value

    */
    int     KP;

    /// index of atom involved in dihedral angle
    /*! AMBER abreviation: LPH, LP, LPPER
        AMBER meaning for LP < 0 is substituted with Type = 1
        index starts from zero and has regular counting order !!
        -1 represents illegal value
    */
    int     LP;

    /// index into parameter arrays DihedralTypes at lambda=1 (initial system)
    /*! AMBER abreviation: ICPH, ICP, ICPPER - first half
        index starts from zero
        -1 represents illegal value
    */
    int     ICP;

    /// index into parameter arrays DihedralTypes at lambda=0 (fully perturbed system)
    /*! AMBER abreviation: ICPPER - second half
        index starts from zero
        -1 represents illegal value
    */
    int     PCP;

    /// type of dihedral angle
    /*! 0 - normal dihedral angle
        1 - improper dihedral angle
       -1 - terminal dihedral angle
    */
    int     Type;

    friend class CAmberDihedralList;
};

//---------------------------------------------------------------------------

#endif
