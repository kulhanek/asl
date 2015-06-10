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
#include <AmberDihedral.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberDihedral::CAmberDihedral(void)
{
    IP = -1;
    JP = -1;
    KP = -1;
    LP = -1;
    ICP = -1;
    PCP = -1;
    Type = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberDihedral::GetIP(void) const
{
    return(IP);
}

//------------------------------------------------------------------------------

void CAmberDihedral::SetIP(int index)
{
    IP = index;
}

//------------------------------------------------------------------------------

int CAmberDihedral::GetJP(void) const
{
    return(JP);
}

//------------------------------------------------------------------------------

void CAmberDihedral::SetJP(int index)
{
    JP = index;
}

//------------------------------------------------------------------------------

int CAmberDihedral::GetKP(void) const
{
    return(KP);
}

//------------------------------------------------------------------------------

void CAmberDihedral::SetKP(int index)
{
    KP = index;
}

//------------------------------------------------------------------------------

int CAmberDihedral::GetLP(void) const
{
    return(LP);
}

//------------------------------------------------------------------------------

void CAmberDihedral::SetLP(int index)
{
    LP = index;
}

//------------------------------------------------------------------------------

int CAmberDihedral::GetICP(void) const
{
    return(ICP);
}

//------------------------------------------------------------------------------

void CAmberDihedral::SetICP(int index)
{
    ICP = index;
}

//------------------------------------------------------------------------------

int CAmberDihedral::GetPCP(void) const
{
    return(PCP);
}

//------------------------------------------------------------------------------

void CAmberDihedral::SetPCP(int index)
{
    PCP = index;
}

//------------------------------------------------------------------------------

int CAmberDihedral::GetType(void) const
{
    return(Type);
}

//------------------------------------------------------------------------------

void CAmberDihedral::SetType(int type)
{
    Type = type;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

