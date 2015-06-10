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
#include <AmberCap.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberCap::CAmberCap(void)
{
    IFCAP = 0;
}

//------------------------------------------------------------------------------

CAmberCap::~CAmberCap(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberCap::InitFields(int iIFCAP)
{
    FreeFields();

    IFCAP = iIFCAP;
    return(true);
}

//------------------------------------------------------------------------------

void CAmberCap::FreeFields(void)
{
    IFCAP = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberCap::HasCap(void) const
{
    return( IFCAP == 1 );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



