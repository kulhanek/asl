#ifndef AmberDihedralListH
#define AmberDihedralListH
/** \ingroup AmberTopology*/
/*! \file AmberDihedralList.hpp */
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
#include <AmberDihedralType.hpp>
#include <AmberDihedral.hpp>

//---------------------------------------------------------------------------

/// dihedral list for topology

class ASL_PACKAGE CAmberDihedralList {
public:
    CAmberDihedralList(void);
    ~CAmberDihedralList(void);

    /// init all fields - old data are destroyed
    void InitFields(int iNPHIH, int iMPHIA, int iNPTRA,
                    int iNPHIA, int iNDPER, int iMDPER);

    /// prepare for new data
    void FreeFields(void);

    /// return number of dihedrals not perturbed
    /*! AMBER abreviation: NPHIH + MPHIA
    */
    int GetNumberOfDihedrals(void) const;

    /// return number of dihedrals containing hydrogen and not perturbed
    /*! AMBER abreviation: NPHIH
    */
    int GetNumberOfDihedralsWithHydrogen(void) const;

    /// return number of dihedrals not containing hydrogen and not perturbed
    /*! AMBER abreviation: MPHIA
    */
    int GetNumberOfDihedralsWithoutHydrogen(void) const;

    /// return number of perturbed dihedrals
    /*! AMBER abreviation: NDPER
    */
    int GetNumberOfPerturbedDihedrals(void) const;

    /// return number of completely perturbed dihedrals
    /*! AMBER abreviation: MDPER
    */
    int GetNumberOfCompletelyPerturbedDihedrals(void) const;

    /// return number of angle types
    /*! AMBER abreviation: NPTRA
    */
    int GetNumberOfDihedralTypes(void) const;

// no longer supported parameters ------------------------------------------
    /// return NPHIA value
    int GetNPHIA(void) const;

// pointers ----------------------------------------------------------------
    /// return pointer to dihedral containing hydrogen
    CAmberDihedral*     GetDihedralWithHydrogen(int index) const;

    /// return pointer to dihedral not containing hydrogen
    CAmberDihedral*     GetDihedralWithoutHydrogen(int index) const;

    /// return pointer to dihedral
    /*!
        method return pointer to any non-perturbed dihedral
        index from
            0 to GetNumberOfDihedralsWithoutHydrogen() - 1 points
            to dihedral not containing hydrogen
        index from
            GetNumberOfDihedralsWithoutHydrogen() to GetNumberOfangles() - 1 points
            to dihedral containing hydrogen
    */
    CAmberDihedral*     GetDihedral(int index) const;

    /// return pointer to perturbed dihedral
    CAmberDihedral*     GetPerturbedDihedral(int index) const;

    /// return pointer to dihedral type
    CAmberDihedralType* GetDihedralType(int index) const;

    /// overload assigment operator
    void operator = (const CAmberDihedralList& src);

    /// remove illegal dihedrals (ICP < 0)
    void RemoveIllegalDihedrals(void);

// section of private data ----------------------------------------------------
private:
    int NPHIH;    //NPHIH  : number of dihedrals containing hydrogen
    int MPHIA;    //MPHIA  : number of dihedrals not containing hydrogen
    int NPTRA;    //NPTRA  : number of unique dihedral types
    int NPHIA;    //NPHIA  : MPHIA + number of constraint dihedrals
    int NDPER;    //NDPER  : number of dihedrals to be perturbed
    int MDPER;    //MDPER  : number of dihedrals with atoms completely in perturbed groups

    CAmberDihedralType* DihedralTypes;
    CAmberDihedral*     DihedralWithoutHydrogens;
    CAmberDihedral*     DihedralWithHydrogens;
    CAmberDihedral*     PerturbedDihedrals;
    bool                SCEEFactorsLoaded;
    bool                SCNBFactorsLoaded;

    bool LoadDihedralPK(FILE* p_file,const char* p_format);
    bool LoadDihedralPN(FILE* p_file,const char* p_format);
    bool LoadDihedralPHASE(FILE* p_file,const char* p_format);
    bool LoadDihedralSCEE(FILE* p_file,const char* p_format);
    bool LoadDihedralSCNB(FILE* p_file,const char* p_format);
    bool LoadDihedralsWithHydrogens(FILE* p_file,const char* p_format);
    bool LoadDihedralsWithoutHydrogens(FILE* p_file,const char* p_format);
    bool LoadPerturbedDihedrals(FILE* p_file,const char* p_format);
    bool LoadPerturbedDihedralTypeIndexes(FILE* p_file,const char* p_format);

    bool SaveDihedralPK(FILE* p_file,const char* p_format);
    bool SaveDihedralPN(FILE* p_file,const char* p_format);
    bool SaveDihedralPHASE(FILE* p_file,const char* p_format);
    bool SaveDihedralSCEE(FILE* p_file,const char* p_format);
    bool SaveDihedralSCNB(FILE* p_file,const char* p_format);
    bool SaveDihedralsWithHydrogens(FILE* p_file,const char* p_format);
    bool SaveDihedralsWithoutHydrogens(FILE* p_file,const char* p_format);
    bool SavePerturbedDihedrals(FILE* p_file,const char* p_format);
    bool SavePerturbedDihedralTypeIndexes(FILE* p_file,const char* p_format);

    friend class CAmberTopology;
    friend class CAmberSubTopology;
};

//---------------------------------------------------------------------------

#endif
