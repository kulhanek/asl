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
#include <AmberAngleType.hpp>
#include <FortranIO.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberAngleType::CAmberAngleType(void)
{
    TK = 0;
    TEQ = 0;
}

//------------------------------------------------------------------------------

double CAmberAngleType::GetTK(void) const
{
    return(TK);
}

//------------------------------------------------------------------------------

void CAmberAngleType::SetTK(double tk)
{
    TK = tk;
}

//------------------------------------------------------------------------------

double CAmberAngleType::GetTEQ(void) const
{
    return(TEQ);
}

//------------------------------------------------------------------------------

void CAmberAngleType::SetTEQ(double teq)
{
    TEQ = teq;
}

//------------------------------------------------------------------------------



