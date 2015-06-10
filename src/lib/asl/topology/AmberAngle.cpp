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
#include <AmberAngle.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberAngle::CAmberAngle(void)
{
    IT = -1;
    JT = -1;
    KT = -1;
    ICT = -1;
    PCT = -1;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberAngle::GetIT(void) const
{
    return(IT);
}

//------------------------------------------------------------------------------

void CAmberAngle::SetIT(int index)
{
    IT = index;
}

//------------------------------------------------------------------------------

int CAmberAngle::GetJT(void) const
{
    return(JT);
}

//------------------------------------------------------------------------------

void CAmberAngle::SetJT(int index)
{
    JT = index;
}


//------------------------------------------------------------------------------

int CAmberAngle::GetKT(void) const
{
    return(KT);
}

//------------------------------------------------------------------------------

void CAmberAngle::SetKT(int index)
{
    KT = index;
}


//------------------------------------------------------------------------------

int CAmberAngle::GetICT(void) const
{
    return(ICT);
}

//------------------------------------------------------------------------------

void CAmberAngle::SetICT(int index)
{
    ICT = index;
}

//------------------------------------------------------------------------------

int CAmberAngle::GetPCT(void) const
{
    return(PCT);
}

//------------------------------------------------------------------------------

void CAmberAngle::SetPCT(int index)
{
    PCT = index;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



