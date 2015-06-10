#ifndef AmberNonBondedListH
#define AmberNonBondedListH
/** \ingroup AmberTopology*/
/*! \file AmberNonBondedList.hpp */
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

#include <stdio.h>
#include <ASLMainHeader.hpp>
#include <AmberAtom.hpp>

//---------------------------------------------------------------------------

/// nonbonded interaction description for topology

class ASL_PACKAGE CAmberNonBondedList {
public:
    CAmberNonBondedList(void);
    ~CAmberNonBondedList(void);

    /// init all fields - previous data are destroyed
    void InitFields(int iNTYPES, int iNEXT, int iNATYP, int iNPHB);

    /// prepare for new data
    void FreeFields(void);

    /// return number of types, e.g. NTYPES
    int  GetNumberOfTypes(void);

    /// return number of excluded atoms, e.g. NEXT
    int GetNumberOfExcludedAtoms(void);

    /// return index of excluded atom
    /*! returned index starts from zero
    */
    int GetNATEX(int index);

    /// set index of excluded atom
    void SetNATEX(int index,int value);

    /// return type of nonbonded interaction
    /*! returned value: sign(ICO[index])
        -1 - 12-10 interaction
        +1 - 12-6 interaction
    */
    int  GetNonBondedType(int index);

    /// return type of nonbonded interaction between two atoms
    /*! returned value:
        -1 - 12-10 interaction
        +1 - 12-6 interaction
    */
    int  GetNonBondedType(CAmberAtom* p_at1,CAmberAtom* p_at2,bool pert=false);

    /// return A parameter of nonbonded interaction between two atoms
    /*! returned value:
        either CN1 or ASOL, type can be determined via GetNonBondedType
    */
    double GetAParam(CAmberAtom* p_at1,CAmberAtom* p_at2,bool pert=false);

    /// return B parameter of nonbonded interaction between two atoms
    /*! returned value:
        either CN2 or BSOL, type can be determined via GetNonBondedType
    */
    double GetBParam(CAmberAtom* p_at1,CAmberAtom* p_at2,bool pert=false);

    /// return icoindex to access nonbonded parameteres
    int GetICOIndex(CAmberAtom* p_at1,CAmberAtom* p_at2,bool pert=false);

    /// return icoindex to access nonbonded parameteres
    int GetICOIndex(int index);

    /// set icoindex to access nonbonded parameteres
    void SetICOIndex(int index,int value);

    /// return A parameter of nonbonded interaction of ico index
    /*! icoindex is non-zero index
        returned value:
        either CN1 (icoindex > 0) or ASOL (icoindex < 0)
    */
    double GetAParam(int icoindex);

    /// return B parameter of nonbonded interaction of ico index
    /*! icoindex is non-zero index
        returned value:
        either CN2 (icoindex > 0) or BSOL (icoindex < 0)
    */
    double GetBParam(int icoindex);

    /// set A parameter of nonbonded interaction of ico index
    /*! icoindex is non-zero index
        method sets
        either CN1 (icoindex > 0) or ASOL (icoindex < 0)
    */
    void SetAParam(double value,int icoindex);

    /// set B parameter of nonbonded interaction of ico index
    /*! icoindex is non-zero index
        method sets
        either CN1 (icoindex > 0) or ASOL (icoindex < 0)
    */
    void SetBParam(double value,int icoindex);

    /// overload assigment operator
    void operator = (const CAmberNonBondedList& src);

// section of private data ----------------------------------------------------
private:
    /// NTYPES : total number of distinct atom types
    int NTYPES;

    /// NEXT   : number of excluded atoms
    int NEXT;

    /// NATYP  : number of atom types in parameter file, see SOLTY below
    int NATYP;

    /// NPHB   : number of distinct 10-12 hydrogen bond pair types
    int NPHB;

    /// ICO    : provides the index to the nonbond parameter arrays CN1, CN2 and ASOL, BSOL
    int*            ICO;
    double*         CN1;
    double*         CN2;
    int*            NATEX;
    double*         ASOL;
    double*         BSOL;
    double*         HBCUT;
    double*         SOLTY;

    bool LoadBasicInfo(FILE* p_file,const char* p_format);
    bool LoadICOs(FILE* p_file,const char* p_format);
    bool LoadSOLTY(FILE* p_file,const char* p_format);
    bool LoadCN1(FILE* p_file,const char* p_format);
    bool LoadCN2(FILE* p_file,const char* p_format);
    bool LoadNATEX(FILE* p_file,const char* p_format);
    bool LoadASOL(FILE* p_file,const char* p_format);
    bool LoadBSOL(FILE* p_file,const char* p_format);
    bool LoadHBCUT(FILE* p_file,const char* p_format);

    bool SaveBasicInfo(FILE* p_file,const char* p_format);
    bool SaveICOs(FILE* p_file,const char* p_format);
    bool SaveSOLTY(FILE* p_file,const char* p_format);
    bool SaveCN1(FILE* p_file,const char* p_format);
    bool SaveCN2(FILE* p_file,const char* p_format);
    bool SaveNATEX(FILE* p_file,const char* p_format);
    bool SaveASOL(FILE* p_file,const char* p_format);
    bool SaveBSOL(FILE* p_file,const char* p_format);
    bool SaveHBCUT(FILE* p_file,const char* p_format);

    friend class CAmberTopology;
};

//---------------------------------------------------------------------------

#endif
