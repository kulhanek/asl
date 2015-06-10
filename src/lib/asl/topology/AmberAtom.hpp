#ifndef AmberAtomH
#define AmberAtomH
/** \ingroup AmberTopology*/
/*! \file AmberAtom.hpp */
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
#include <set>

// -----------------------------------------------------------------------------

class CAmberResidue;

// -----------------------------------------------------------------------------

/// atom description

class ASL_PACKAGE CAmberAtom {
public:
    CAmberAtom(void);

    /// return number of excluded atoms with this atom
    int    GetNUMEX(void) const;

    /// set number of excluded atoms with thos atom
    void   SetNUMEX(int iNUMEX);

    /// return atom mass
    /*! unit of mass : g/mol
    */
    double GetMass(void) const;

    /// set atom mass
    /*! unit of mass : g/mol
    */
    void SetMass(double mass);

    /// return name of atom
    const char*  GetName(bool pert=false) const;

    /// set name of atom
    void  SetName(const char* p_name,bool pert=false);

    /// return tree name of atom
    const char*  GetTreeName(void) const;

    /// set tree name of atom
    void  SetTreeName(const char* p_tname);

    /// return name in PDB format
    const CSmallString  GetPDBName(bool pert=false) const;

    /// return type of atom
    const char*  GetType(bool pert=false) const;

    /// set type of atom
    void  SetType(const char* p_type,bool pert=false);

    /// return radius of atom
    double GetRadius(void) const;

    /// set radius of atom
    void SetRadius(double radius);

    /// return screen value of atom
    double GetScreenValue(void) const;

    /// set screen value of atom
    void SetScreenValue(double screen);

    /// return charge of atom
    /*! unit of charge : internal amber charge unit
    */
    double GetCharge(bool pert=false) const;

    /// set charge
    /*! unit of charge : internal amber charge unit
    */
    void SetCharge(double charge,bool pert=false);

    /// return charge of atom
    /*! unit of charge : units of electron charge
        this value is recalculated from AMBER internal value !!!
        e.g. it is divided by 18.2223
    */
    double GetStandardCharge(bool pert=false) const;

    /// return charge of atom
    /*! unit of charge : units of electron charge
        this value is recalculated to AMBER internal value !!!
        e.g. it is multiplied by 18.2223
    */
    void SetStandardCharge(double charge,bool pert=false);

    /// return IAC index for nonbonded interaction
    int    GetIAC(bool pert=false) const;

    /// set IAC index for nonbonded interaction
    void    SetIAC(int iac_index,bool pert=false);

    /// return atom polarizability
    double GetPol(bool pert=false) const;

    /// set atom polarizability
    void SetPol(double pol,bool pert=false);

    /// return atomic number
    int GetAtomicNumber(void) const;

    /// set atomic number
    void SetAtomicNumber(int z);

    /// return true if atom is perturbed
    bool IsPerturbed(void) const;

    /// return pointer to residue in which atom belongs
    CAmberResidue* GetResidue(void) const;

    /// return index of molecule in which atom belongs
    /*! this information is initialized by CAmberTopology:InitMoleculeIndexes()
        -1 represents non-initialized state, index does not have to be continuous
    */
    int GetMoleculeIndex(void) const;

    /// return index of atom
    int GetAtomIndex(void) const;

    /// guess atom number (z)
    int GuessZ(void) const;

    /// get number of neighbour atoms
    int GetNumberOfNeighbourAtoms(void);

    /// get neighbour atom index
    int GetNeighbourAtomIndex(int index);

// section of private data ----------------------------------------------------
private:
    // common properties -------------------------------------------------------

    /// NUMEX  : total number of excluded atoms for atom "i".  See NATEX below.
    int     NUMEX;

    /// ITREE  :  the list of tree joining information
    /*! classified into five types
        M -- main chain,
        S -- side chain,
        B -- branch point,
        3 -- branch into three chains,
        E -- end of the chain
    */
    char    ITREE[5];

    /// JOIN : tree joining information, potentially used in ancient analysis programs.
    int     JOIN;

    /// IROTAT : no longer used
    int     IROTAT;

    /// AMASS  : the atom masses
    /*! unit of mass : g/mol
    */
    double  AMASS;

    /// RADIUS - only if ver >= 7
    double  RADIUS;

    /// SCREEN - only if ver >= 7
    double  SCREEN;

    // properties at lambda=1 (initial state) ----------------------------------

    /// IGRAPH : the user atoms names - zero terminated string
    char    IGRAPH[5];

    /// ISYMBL : the AMBER atom types for each atom - zero terminated string
    char    ISYMBL[5];

    /// CHRG   : the atom charges
    /*! unit of charge : units of electron charge
        WARNING: this value is recalculated from AMBER internal value !!!
        e.g. it is divided by 18.2223
    */
    double  CHRG;

    /// IAC    : index for the atom types involved in Lennard Jones (6-12) interactions.
    int     IAC;

    /// ATPOL  : atomic polarizabilities
    double  ATPOL;

    // properties at lambda=0 (fully perturbed state) --------------------------

    /// IAPER  : IAPER = 1 if the atom is being perturbed
    int     IAPER;

    /// IGRPER : atomic names at lambda=0, zero terminated string
    char    IGRPER[5];

    /// ISMPER : atomic symbols at lambda=0, zero terminated string
    char    ISMPER[5];

    /// CGPER  : atomic charges at lambda=0
    /*! unit of charge : units of electron charge
        WARNING: this value is recalculated from AMBER internal value !!!
        e.g. it is divided by 18.2223
    */
    double  CGPER;

    /// IACPER : index for the atom types involved in Lennard Jones interactions at lambda=0
    int     IACPER;

    /// ATPOL1 : atomic polarizabilities at lambda = 0
    double  ATPOL1;

    /// no longer used value
    double  ALMPER;

    /// atomic number
    int     ATOMIC_NUMBER;

    // additional information --------------------------------------------------
    /// pointer to residue in which atom belongs
    CAmberResidue*  Residue;

    /// molecule index
    int             MoleculeIndex;

    /// atom index
    int             AtomIndex;

    /// list of neighbour atom indexes
    std::set<int>   NeighbourIndexes;

    friend class CAmberAtomList;
    friend class CAmberResidueList;
    friend class CAmberTopology;
};

// -----------------------------------------------------------------------------

#endif
