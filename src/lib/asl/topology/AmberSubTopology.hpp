#ifndef AmberSubTopologyH
#define AmberSubTopologyH
/** \ingroup AmberTopology*/
/*! \file AmberSubTopology.hpp */
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
#include <AmberTopology.hpp>

//---------------------------------------------------------------------------

class CAmberTopology;
class CAmberMaskAtoms;

//---------------------------------------------------------------------------

/// enable to create sub topology

class ASL_PACKAGE CAmberSubTopology : public CAmberTopology {
public:
    CAmberSubTopology(void);
    ~CAmberSubTopology(void);

    /// init subtopology
    bool InitSubTopology(CAmberMaskAtoms* p_mask,
                                      bool copy_box = false,
                                      bool ignore_errors = false,
                                      bool verbose = false);

    /// set stream for report
    void SetReportFile(FILE* p_fout);

// section of private data ----------------------------------------------------
private:
    CAmberTopology*     OldTopology;
    CAmberMaskAtoms*    Mask;

    FILE*               ReportFile;     // file for report
    bool                VerboseMode;    // print informative messages
    bool                IgnoreErrors;   // ignore errors
    bool                CopyBox;        // copy box info

    int*                AtomMapper;     // atom mapper for atom reindexing
    bool                CrossBonds;     // were cross-bonds detected?

    bool PrepareNewTopology(void);
    bool PrepareAtoms(void);
    bool PrepareNonbondedList(void);
    bool PrepareResidues(void);
    bool PrepareBonds(void);
    bool PrepareAngles(void);
    bool PrepareDihedrals(void);
    bool CopyTopologyBox(void);
};

//---------------------------------------------------------------------------

#endif
