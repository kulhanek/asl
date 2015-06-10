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
#include <AmberDihedralType.hpp>
#include <FortranIO.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberDihedralType::CAmberDihedralType(void)
{
    PK = 0;
    PN = 0;
    PHASE = 0;
    // default AMBER factors
    SCEE_SCALE = 1.2;
    SCNB_SCALE = 2.0;
}

//------------------------------------------------------------------------------

double CAmberDihedralType::GetPK(void) const
{
    return(PK);
}

//------------------------------------------------------------------------------

void CAmberDihedralType::SetPK(double pk)
{
    PK = pk;
}

//------------------------------------------------------------------------------

double CAmberDihedralType::GetPN(void) const
{
    return(PN);
}

//------------------------------------------------------------------------------

void CAmberDihedralType::SetPN(double pn)
{
    PN = pn;
}

//------------------------------------------------------------------------------

double CAmberDihedralType::GetPHASE(void) const
{
    return(PHASE);
}

//------------------------------------------------------------------------------

void CAmberDihedralType::SetPHASE(double phase)
{
    PHASE = phase;
}

//------------------------------------------------------------------------------

double CAmberDihedralType::GetSCEE(void) const
{
    return(SCEE_SCALE);
}

//------------------------------------------------------------------------------

void CAmberDihedralType::SetSCEE(double scee)
{
    SCEE_SCALE = scee;
}

//------------------------------------------------------------------------------

double CAmberDihedralType::GetSCNB(void) const
{
    return(SCNB_SCALE);
}

//------------------------------------------------------------------------------

void CAmberDihedralType::SetSCNB(double scnb)
{
    SCNB_SCALE = scnb;
}

//------------------------------------------------------------------------------



