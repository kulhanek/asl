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
#include <AmberResidue.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberResidue::CAmberResidue(void)
{
    Topology = NULL;
    memset(LABRES,0,5);
    memset(PERRES,0,5);
    IPRES=0;
    Index = -1;
    NumOfAtoms = 0;
    NumOfBondsWithHydrogen = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

const char* CAmberResidue::GetName(bool pert) const
{
    if( pert == false ) {
        return(LABRES);
    } else {
        return(PERRES);
    }
}

//------------------------------------------------------------------------------

void  CAmberResidue::SetName(const char* p_name,bool pert)
{
    if( pert == false ) {
        strncpy(LABRES,p_name,4);
    } else {
        strncpy(PERRES,p_name,4);
    }
}

//------------------------------------------------------------------------------

int CAmberResidue::GetNumberOfAtoms(void) const
{
    return(NumOfAtoms);
}

//------------------------------------------------------------------------------

int CAmberResidue::GetFirstAtomIndex(void) const
{
    return(IPRES-1);
}

//------------------------------------------------------------------------------

void  CAmberResidue::SetFirstAtomIndex(int index)
{
    IPRES = index + 1;
}

//------------------------------------------------------------------------------

int CAmberResidue::GetIndex(void) const
{
    return(Index);
}

//------------------------------------------------------------------------------

int CAmberResidue::GetNumberOfBondsWithHydrogen(void) const
{
    return(NumOfBondsWithHydrogen);
}

//------------------------------------------------------------------------------

CAmberTopology* CAmberResidue::GetTopology(void) const
{
    return(Topology);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




