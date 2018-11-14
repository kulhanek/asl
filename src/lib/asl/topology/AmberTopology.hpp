#ifndef topologyH
#define topologyH
/** \ingroup AmberTopology*/
/*! \file AmberTopology.hpp */
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
#include <AmberAtomList.hpp>
#include <AmberResidueList.hpp>
#include <AmberBondList.hpp>
#include <AmberAngleList.hpp>
#include <AmberDihedralList.hpp>
#include <AmberNonBondedList.hpp>
#include <AmberBox.hpp>
#include <AmberCap.hpp>

//---------------------------------------------------------------------------

/// amber topology versions

enum EAmberVersion {
    AMBER_VERSION_NONE,
    AMBER_VERSION_6,
    AMBER_VERSION_7
};

//---------------------------------------------------------------------------

/// topology description

class ASL_PACKAGE CAmberTopology {
public:
    CAmberTopology(void);
    ~CAmberTopology(void);

    /// load topology, if allow_stdin==true then file_name == '-' means stdin
    bool Load(const CSmallString& file_name,bool allow_stdin=false);

    /// load topology from a file stream
    bool Load(FILE* p_fin);

    /// load fake topology, if allow_stdin==true then file_name == '-' means stdin
    bool LoadFakeTopologyFromPDB(const CSmallString& file_name,bool mangle_names,bool allow_stdin=false);

    /// load fake topology from a file stream
    bool LoadFakeTopologyFromPDB(FILE* p_fin,bool mangle_names);

    /// load fake topology, if allow_stdin==true then file_name == '-' means stdin
    bool LoadFakeTopologyFromG96(const CSmallString& file_name,bool allow_stdin=false);

    /// load fake topology from a file stream
    bool LoadFakeTopologyFromG96(FILE* p_fin);

    /// is fake topology?
    bool IsFake(void);

    /// save topology - if version is not specified, version from load is used
    bool Save(const CSmallString& file_name,
                           EAmberVersion version=AMBER_VERSION_NONE,
                           bool allow_stdout=false);

    /// save topology - if version is not specified, version from load is used
    bool Save(FILE* p_fout,EAmberVersion version=AMBER_VERSION_NONE);

    /// prepare for new data
    void Clean(void);

    /// print info about topology
    void PrintInfo(bool short_info=false,FILE* p_out=NULL);

    /// print info about atoms
    void PrintAtoms(FILE* p_out=NULL);

    /// return version of currently loaded topology
    EAmberVersion GetVersion(void);

    /// return title of topology
    CSmallString GetTitle(void);

    /// set title of topology
    void SetTitle(const CSmallString& title);

    /// init GetNumberOfBondsWithHydrogen for CAmberResidue
    bool InitResidueNumOfBondsWithHydrogen(void);

    /// init id of molecules for individual atoms, it is not influenced by informations in BoxInfo
    bool InitMoleculeIndexes(void);

    /// build list of neighbour atoms
    void BuidListOfNeighbourAtoms(void);

    /// overload assigment operator
    void operator = (const CAmberTopology& src);
    
    /// get total mass
    double GetTotalMass(void);

// section o public data ------------------------------------------------------
public:
    CAmberAtomList      AtomList;
    CAmberResidueList   ResidueList;
    CAmberBondList      BondList;
    CAmberAngleList     AngleList;
    CAmberDihedralList  DihedralList;
    CAmberNonBondedList NonBondedList;
    CAmberBox           BoxInfo;
    CAmberCap           CapInfo;
    bool                FakeTopology;

// section of private data ----------------------------------------------------
private:
    /// currently loaded version of topology
    EAmberVersion   Version;

    /// name of topology
    CSmallString    Name;

    /// ITITL  : title
    CSmallString    ITITL;

    /// NHPARM : currently not used
    int NHPARM;

    /// NPARM  : currently not used
    int NPARM;
    
    // total mass in g/mol
    double      TotalMass;

    // local copy of formats
    CSmallString fTITLE;
    CSmallString fPOINTERS;
    CSmallString fATOM_NAME;
    CSmallString fCHARGE;
    CSmallString fMASS;
    CSmallString fATOM_TYPE_INDEX;
    CSmallString fNUMBER_EXCLUDED_ATOMS;
    CSmallString fNONBONDED_PARM_INDEX;
    CSmallString fRESIDUE_LABEL;
    CSmallString fRESIDUE_POINTER;
    CSmallString fBOND_FORCE_CONSTANT;
    CSmallString fBOND_EQUIL_VALUE;
    CSmallString fANGLE_FORCE_CONSTANT;
    CSmallString fANGLE_EQUIL_VALUE;
    CSmallString fDIHEDRAL_FORCE_CONSTANT;
    CSmallString fDIHEDRAL_PERIODICITY;
    CSmallString fDIHEDRAL_PHASE;
    CSmallString fSOLTY;
    CSmallString fLENNARD_JONES_ACOEF;
    CSmallString fLENNARD_JONES_BCOEF;
    CSmallString fBONDS_INC_HYDROGEN;
    CSmallString fBONDS_WITHOUT_HYDROGEN;
    CSmallString fANGLES_INC_HYDROGEN;
    CSmallString fANGLES_WITHOUT_HYDROGEN;
    CSmallString fDIHEDRALS_INC_HYDROGEN;
    CSmallString fDIHEDRALS_WITHOUT_HYDROGEN;
    CSmallString fEXCLUDED_ATOMS_LIST;
    CSmallString fHBOND_ACOEF;
    CSmallString fHBOND_BCOEF;
    CSmallString fHBCUT;
    CSmallString fAMBER_ATOM_TYPE;
    CSmallString fTREE_CHAIN_CLASSIFICATION;
    CSmallString fJOIN_ARRAY;
    CSmallString fIROTAT;
    CSmallString fSOLVENT_POINTERS;
    CSmallString fATOMS_PER_MOLECULE;
    CSmallString fBOX_DIMENSIONS;
    CSmallString fRADIUS_SET;
    CSmallString fRADII;
    CSmallString fSCREEN;
    CSmallString fPERT_BOND_ATOMS;
    CSmallString fPERT_BOND_PARAMS;
    CSmallString fPERT_ANGLE_ATOMS;
    CSmallString fPERT_ANGLE_PARAMS;
    CSmallString fPERT_DIHEDRAL_ATOMS;
    CSmallString fPERT_DIHEDRAL_PARAMS;
    CSmallString fPERT_RESIDUE_NAME;
    CSmallString fPERT_ATOM_NAME;
    CSmallString fPERT_ATOM_SYMBOL;
    CSmallString fALMPER;
    CSmallString fIAPER;
    CSmallString fPERT_ATOM_TYPE_INDEX;
    CSmallString fPERT_CHARGE;
    CSmallString fSCEE_SCALE_FACTOR;
    CSmallString fSCNB_SCALE_FACTOR;
    CSmallString fATOMIC_NUMBER;
    CSmallString fIPOL;
    CSmallString fPOL;

    bool LoadAmber6(FILE* p_top);
    bool LoadAmber7(FILE* p_top);
    bool SaveAmber6(FILE* p_top);
    bool SaveAmber7(FILE* p_top);

    bool LoadBasicInfo(FILE* p_top,const char* p_format,EAmberVersion version);
    bool SaveBasicInfo(FILE* p_top,const char* p_format,EAmberVersion version);
    bool SaveSectionHeader(FILE* p_top,const char* p_section_name,
                                        const char* p_section_format);

    void SetDefaultAmber7Formats(void);

private:
    // disable copy constructor
    CAmberTopology(const CAmberTopology& copy);
};

//---------------------------------------------------------------------------

#endif
