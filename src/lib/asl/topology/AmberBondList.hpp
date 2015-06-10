#ifndef AmberBondListH
#define AmberBondListH
/** \ingroup AmberTopology*/
/*! \file AmberBondList.hpp */
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
#include <AmberBondType.hpp>
#include <AmberBond.hpp>

//---------------------------------------------------------------------------

/// list of bonds for topology
class ASL_PACKAGE CAmberBondList {
public:
    CAmberBondList(void);
    ~CAmberBondList(void);

// public methods ----------------------------------------------------------
    /// init all fields - old data are destroyed
    void InitFields(int iNBONH,int iMBONA,int iNUMBND,int iNBONA,
                    int iNBPER,int iMBPER);

    /// prepare for new data
    void FreeFields(void);

    /// return number of bonds not perturbed
    /*! AMBER abreviation: MBONA + NBONH
    */
    int GetNumberOfBonds(void) const;

    /// return number of bonds containing hydrogen and not perturbed
    /*! AMBER abreviation: NBONH
    */
    int GetNumberOfBondsWithHydrogen(void) const;

    /// return number of bonds not containing hydrogen and not perturbed
    /*! AMBER abreviation: MBONA
    */
    int GetNumberOfBondsWithoutHydrogen(void) const;

    /// return number of perturbed bonds
    /*! AMBER abreviation: NBPER
    */
    int GetNumberOfPerturbedBonds(void) const;

    /// return number of completely perturbed bonds
    /*! AMBER abreviation: MBPER
    */
    int GetNumberOfCompletelyPerturbedBonds(void) const;

    /// return number of bond types
    /*! AMBER abreviation: NUMBND
    */
    int GetNumberOfBondTypes(void) const;

// no longer supported parameters ------------------------------------------
    /// return NBONA value
    int GetNBONA(void) const;

// pointers ----------------------------------------------------------------
    /// return pointer to bond containing hydrogen
    CAmberBond*     GetBondWithHydrogen(int index) const;

    /// return pointer to bond not containing hydrogen
    CAmberBond*     GetBondWithoutHydrogen(int index) const;

    /// return pointer to bond
    /*!
        method return pointer to any non-perturbed bond
        index from 0 to GetNumberOfBondsWithoutHydrogen() - 1 points to bond not containing hydrogen
        index from GetNumberOfBondsWithoutHydrogen() to GetNumberOfBonds() - 1 points to bond containing hydrogen
    */
    CAmberBond*     GetBond(int index) const;

    /// return pointer to perturbed bond
    CAmberBond*     GetPerturbedBond(int index) const;

    /// return pointer to bond type
    CAmberBondType* GetBondType(int index) const;

    /// overload assigment operator
    void operator = (const CAmberBondList& src);

// section of private data ----------------------------------------------------
private:
    int     NBONH;    //NBONH  : number of bonds containing hydrogen
    int     MBONA;    //MBONA  : number of bonds not containing hydrogen
    int     NUMBND;   //NUMBND : number of unique bond types
    int     NBONA;    //NBONA  : MBONA + number of constraint bonds
    int     NBPER;    //NBPER  : number of bonds to be perturbed
    int     MBPER;    //MBPER  : number of bonds with atoms completely in perturbed group

    CAmberBondType* BondTypes;
    CAmberBond*     BondsWithoutHydrogens;
    CAmberBond*     BondsWithHydrogens;
    CAmberBond*     PerturbedBonds;

    bool LoadBondRK(FILE* p_file,const char* p_format);
    bool LoadBondREQ(FILE* p_file,const char* p_format);
    bool LoadBondsWithHydrogens(FILE* p_file,const char* p_format);
    bool LoadBondsWithoutHydrogens(FILE* p_file,const char* p_format);
    bool LoadPerturbedBonds(FILE* p_file,const char* p_format);
    bool LoadPerturbedBondTypeIndexes(FILE* p_file,const char* p_format);

    bool SaveBondRK(FILE* p_file,const char* p_format);
    bool SaveBondREQ(FILE* p_file,const char* p_format);
    bool SaveBondsWithHydrogens(FILE* p_file,const char* p_format);
    bool SaveBondsWithoutHydrogens(FILE* p_file,const char* p_format);
    bool SavePerturbedBonds(FILE* p_file,const char* p_format);
    bool SavePerturbedBondTypeIndexes(FILE* p_file,const char* p_format);

    friend class CAmberTopology;
};

//---------------------------------------------------------------------------
#endif
