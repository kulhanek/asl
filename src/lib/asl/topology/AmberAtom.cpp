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

#include <string.h>
#include <stdlib.h>
#include <AmberAtom.hpp>

#include <PeriodicTable.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberAtom::CAmberAtom(void)
{
    memset(IGRAPH,0,5);
    CHRG = 0;
    AMASS = 0;
    IAC = 0;
    NUMEX = 0;
    memset(ISYMBL,0,5);
    memset(ITREE,0,5);
    JOIN = 0;
    IROTAT = 0;
    memset(IGRPER,0,5);
    memset(ISMPER,0,5);
    ALMPER = 0;
    IAPER = 0;
    IACPER = 0;
    CGPER = 0;
    ATPOL = 0;
    ATPOL1 = 0;
    RADIUS = 0;
    SCREEN = 0;
    ATOMIC_NUMBER = 0;
    Residue = NULL;
    MoleculeIndex = -1;
    AtomIndex = -1;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberAtom::GetNUMEX(void) const
{
    return(NUMEX);
}

//------------------------------------------------------------------------------

void CAmberAtom::SetNUMEX(int iNUMEX)
{
    NUMEX = iNUMEX;
}

//------------------------------------------------------------------------------

double CAmberAtom::GetMass(void) const
{
    return(AMASS);
}

//------------------------------------------------------------------------------

void CAmberAtom::SetMass(double mass)
{
    AMASS = mass;
}

//------------------------------------------------------------------------------

const char* CAmberAtom::GetName(bool pert) const
{
    if( pert == false ) {
        return(IGRAPH);
    } else {
        return(IGRPER);
    }
}

//------------------------------------------------------------------------------

void  CAmberAtom::SetName(const char* p_name,bool pert)
{
    if( pert == false ) {
        strncpy(IGRAPH,p_name,4);
    } else {
        strncpy(IGRPER,p_name,4);
    }
}

//------------------------------------------------------------------------------

const char* CAmberAtom::GetTreeName(void) const
{
    return(ITREE);
}

//------------------------------------------------------------------------------

void  CAmberAtom::SetTreeName(const char* p_tname)
{
    strncpy(ITREE,p_tname,4);
}

//------------------------------------------------------------------------------

const CSmallString  CAmberAtom::GetPDBName(bool pert) const
{
    char buffer[5];
    buffer[4]='\0';
    buffer[0] = GetName(pert)[3];
    buffer[1] = GetName(pert)[0];
    buffer[2] = GetName(pert)[1];
    buffer[3] = GetName(pert)[2];
    return(CSmallString(buffer));
}

//------------------------------------------------------------------------------

const char* CAmberAtom::GetType(bool pert) const
{
    if( pert == false ) {
        return(ISYMBL);
    } else {
        return(ISMPER);
    }
}

//------------------------------------------------------------------------------

void CAmberAtom::SetType(const char* p_type,bool pert)
{
    if( pert == false ) {
        strncpy(ISYMBL,p_type,4);
    } else {
        strncpy(ISMPER,p_type,4);
    }
}

//------------------------------------------------------------------------------

double CAmberAtom::GetRadius(void) const
{
    return(RADIUS);
}

//------------------------------------------------------------------------------

void CAmberAtom::SetRadius(double radius)
{
    RADIUS = radius;
}

//------------------------------------------------------------------------------

double CAmberAtom::GetScreenValue(void) const
{
    return(SCREEN);
}

//------------------------------------------------------------------------------

void CAmberAtom::SetScreenValue(double screen)
{
    SCREEN = screen;
}

//------------------------------------------------------------------------------

double CAmberAtom::GetCharge(bool pert) const
{
    if( pert == false ) {
        return(CHRG);
    } else {
        return(CGPER);
    }
}

//------------------------------------------------------------------------------

void CAmberAtom::SetCharge(double charge,bool pert)
{
    if( pert == false ) {
        CHRG = charge;
    } else {
        CGPER = charge;
    }
}

//------------------------------------------------------------------------------

double CAmberAtom::GetStandardCharge(bool pert) const
{
    if( pert == false ) {
        return(CHRG/18.2223);
    } else {
        return(CGPER/18.2223);
    }
}

//------------------------------------------------------------------------------

void CAmberAtom::SetStandardCharge(double charge,bool pert)
{
    if( pert == false ) {
        CHRG = charge*18.2223;
    } else {
        CGPER = charge*18.2223;
    }
}

//------------------------------------------------------------------------------

int CAmberAtom::GetIAC(bool pert) const
{
    if( pert == false ) {
        return(IAC);
    } else {
        return(IACPER);
    }
}

//------------------------------------------------------------------------------

void CAmberAtom::SetIAC(int iac_index,bool pert)
{
    if( pert == false ) {
        IAC = iac_index;
    } else {
        IACPER = iac_index;
    }
}

//------------------------------------------------------------------------------

double CAmberAtom::GetPol(bool pert) const
{
    if( pert == false ) {
        return(ATPOL);
    } else {
        return(ATPOL1);
    }
}

//------------------------------------------------------------------------------

void CAmberAtom::SetPol(double pol,bool pert)
{
    if( pert == false ) {
        ATPOL = pol;
    } else {
        ATPOL1 = pol;
    }
}

//------------------------------------------------------------------------------

int CAmberAtom::GetAtomicNumber(void) const
{
    return(ATOMIC_NUMBER);
}

//------------------------------------------------------------------------------

void CAmberAtom::SetAtomicNumber(int z)
{
    ATOMIC_NUMBER = z;
}

//------------------------------------------------------------------------------

bool CAmberAtom::IsPerturbed(void) const
{
    return( IAPER == 1 );
}

//------------------------------------------------------------------------------

CAmberResidue* CAmberAtom::GetResidue(void) const
{
    return(Residue);
}

//------------------------------------------------------------------------------

int CAmberAtom::GetMoleculeIndex(void) const
{
    return(MoleculeIndex);
}

//------------------------------------------------------------------------------

int CAmberAtom::GetAtomIndex(void) const
{
    return(AtomIndex);
}

//------------------------------------------------------------------------------

int CAmberAtom::GuessZ(void) const
{
    if( ATOMIC_NUMBER != 0 ) return(ATOMIC_NUMBER);
    int z = PeriodicTable.SearchZByMass(this->GetMass());
    return (z);
}

//------------------------------------------------------------------------------

int CAmberAtom::GetNumberOfNeighbourAtoms(void)
{
    return(NeighbourIndexes.size());
}

//------------------------------------------------------------------------------

int CAmberAtom::GetNeighbourAtomIndex(int index)
{
    std::set<int>::iterator it = NeighbourIndexes.begin();
    for(int i=0; i < index; i++) it++;
    if( it != NeighbourIndexes.end() ){
        return(*it);
    } else {
        return(-1);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




