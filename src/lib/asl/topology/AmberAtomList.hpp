#ifndef AmberAtomListH
#define AmberAtomListH
/** \ingroup AmberTopology*/
/*! \file AmberAtomList.hpp */
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
#include <SmallString.hpp>
#include <AmberAtom.hpp>

//---------------------------------------------------------------------------

/// atom list for topology

class ASL_PACKAGE CAmberAtomList {
public:
    CAmberAtomList(void);
    ~CAmberAtomList(void);

    /// init all fields - old data are destroyed
    void InitFields(int iNATOM,int iIFPERT,int iIFPOL);

    /// prepare for new data
    void FreeFields(void);

    /// return number of all atoms
    int  GetNumberOfAtoms(void) const;

    /// return pointer to atom
    CAmberAtom* GetAtom(int index) const;

    /// return true if the system contains perturbation information
    bool HasPertInfo(void) const;

    /// return true if the system contains polarization information
    bool HasPolInfo(void) const;

    /// get name of radius set
    const CSmallString& GetRadiusSet(void) const;

    /// set name of radisu set
    void SetRadiusSet(const CSmallString& set_name);

    /// overload assigment operator
    void operator = (const CAmberAtomList& src);

// section of private data ----------------------------------------------------
private:
    /// NATOM  : total number of atoms
    int         NATOM;

    /// IFPERT : set to 1 if perturbation info is to be read in
    int         IFPERT;

    /// IFPOL : set to 1 if polarization info is to be read in
    int         IFPOL;

    // options
    bool            AtomicNumberLoaded;
    CSmallString    RadiusSet;

    CAmberAtom* Atoms;

    bool LoadAtomNames(FILE* p_file,const char* p_format);
    bool LoadAtomCharges(FILE* p_file,const char* p_format);
    bool LoadAtomAtomicNumbers(FILE* p_file,const char* p_format);
    bool LoadAtomMasses(FILE* p_file,const char* p_format);
    bool LoadAtomIACs(FILE* p_file,const char* p_format);
    bool LoadAtomNUMEXs(FILE* p_file,const char* p_format);
    bool LoadAtomIPol(FILE* p_file,const char* p_format);
    bool LoadAtomPol(FILE* p_file,const char* p_format);

    bool SaveAtomNames(FILE* p_file,const char* p_format);
    bool SaveAtomCharges(FILE* p_file,const char* p_format);
    bool SaveAtomAtomicNumbers(FILE* p_file,const char* p_format);
    bool SaveAtomMasses(FILE* p_file,const char* p_format);
    bool SaveAtomIACs(FILE* p_file,const char* p_format);
    bool SaveAtomNUMEXs(FILE* p_file,const char* p_format);
    bool SaveAtomPol(FILE* p_file,const char* p_format);

    bool LoadAtomISYMBL(FILE* p_file,const char* p_format);
    bool LoadAtomITREE(FILE* p_file,const char* p_format);
    bool LoadAtomJOIN(FILE* p_file,const char* p_format);
    bool LoadAtomIROTAT(FILE* p_file,const char* p_format);

    bool SaveAtomISYMBL(FILE* p_file,const char* p_format);
    bool SaveAtomITREE(FILE* p_file,const char* p_format);
    bool SaveAtomJOIN(FILE* p_file,const char* p_format);
    bool SaveAtomIROTAT(FILE* p_file,const char* p_format);

    bool LoadPertAtomNames(FILE* p_file,const char* p_format);
    bool LoadPertAtomISYMBL(FILE* p_file,const char* p_format);
    bool LoadPertAtomCharges(FILE* p_file,const char* p_format);
    bool LoadPertAtomPertFlag(FILE* p_file,const char* p_format);
    bool LoadPertAtomIAC(FILE* p_file,const char* p_format);
    bool LoadPertAtomPol(FILE* p_file,const char* p_format);
    bool LoadPertAtomALMPER(FILE* p_file,const char* p_format);

    bool SavePertAtomNames(FILE* p_file,const char* p_format);
    bool SavePertAtomISYMBL(FILE* p_file,const char* p_format);
    bool SavePertAtomCharges(FILE* p_file,const char* p_format);
    bool SavePertAtomPertFlag(FILE* p_file,const char* p_format);
    bool SavePertAtomIAC(FILE* p_file,const char* p_format);
    bool SavePertAtomPol(FILE* p_file,const char* p_format);
    bool SavePertAtomALMPER(FILE* p_file,const char* p_format);

    bool LoadAtomRadiusSet(FILE* p_file,const char* p_format);
    bool LoadAtomRadii(FILE* p_file,const char* p_format);
    bool LoadAtomScreen(FILE* p_file,const char* p_format);
    bool SaveAtomRadiusSet(FILE* p_file,const char* p_format);
    bool SaveAtomRadii(FILE* p_file,const char* p_format);
    bool SaveAtomScreen(FILE* p_file,const char* p_format);

    friend class CAmberTopology;
    friend class CAmberSubTopology;
};

//---------------------------------------------------------------------------

#endif
