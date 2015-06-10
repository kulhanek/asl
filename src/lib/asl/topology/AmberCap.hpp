#ifndef AmberCapH
#define AmberCapH
/** \ingroup AmberTopology*/
/*! \file AmberCap.hpp */
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

class CAmberResidue;

//---------------------------------------------------------------------------

/// CAP description

class ASL_PACKAGE CAmberCap {
public:
    CAmberCap(void);
    ~CAmberCap(void);

    /// init all fields - old data are destroyed
    bool InitFields(int iIFCAP);

    /// prepare for new data
    void FreeFields(void);

    /// return true if topology contains info about CAP
    bool HasCap(void) const;

// section of private data ----------------------------------------------------
private:
    /// IFCAP : set to 1 if the CAP option from edit was specified
    int     IFCAP;

    friend class CAmberTopology;
};

//---------------------------------------------------------------------------

#endif
