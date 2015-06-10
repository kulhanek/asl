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
#include <AmberBondType.hpp>
#include <FortranIO.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberBondType::CAmberBondType(void)
{
    RK = 0;
    REQ = 0;
}

//------------------------------------------------------------------------------

double CAmberBondType::GetRK(void) const
{
    return(RK);
}

//------------------------------------------------------------------------------

void CAmberBondType::SetRK(double rk)
{
    RK = rk;
}

//------------------------------------------------------------------------------

double CAmberBondType::GetREQ(void) const
{
    return(REQ);
}

//------------------------------------------------------------------------------

void CAmberBondType::SetREQ(double req)
{
    REQ = req;
}

// -----------------------------------------------------------------------------

