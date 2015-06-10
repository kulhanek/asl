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
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include <AmberMaskAtoms.hpp>
#include <AmberTopology.hpp>
#include <AmberRestart.hpp>
#include <FortranIO.hpp>
#include <ErrorSystem.hpp>
#include <AmberMaskASelection.hpp>

#include "maskparser/AmberMaskParser.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberMaskAtoms::CAmberMaskAtoms(void)
{
    Topology = NULL;
    Coordinates = NULL;
    Selection = NULL;
}

//------------------------------------------------------------------------------

CAmberMaskAtoms::~CAmberMaskAtoms(void)
{
    if( Selection != NULL ) delete Selection;
    Selection = NULL;
    Topology = NULL;
    Coordinates = NULL;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberMaskAtoms::AssignTopology(CAmberTopology* p_top)
{
    if( Selection != NULL ) delete Selection;
    SelectedAtoms.clear();
    Topology = p_top;
    Selection = NULL;
    Coordinates = NULL;
    Mask = NULL;
    return(true);
}

//------------------------------------------------------------------------------

CAmberTopology* CAmberMaskAtoms::GetTopology(void) const
{
    return(Topology);
}

//------------------------------------------------------------------------------

bool CAmberMaskAtoms::AssignCoordinates(CAmberRestart* p_crd)
{
    Coordinates = p_crd;
    return(true);
}

//------------------------------------------------------------------------------

CAmberRestart* CAmberMaskAtoms::GetCoordinates(void) const
{
    return(Coordinates);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberMaskAtoms::SelectAllAtoms(void)
{
    return(SetMask("@*"));
}

//------------------------------------------------------------------------------

void CAmberMaskAtoms::Reset(void)
{
    if( Selection != NULL ) delete Selection;
    SelectedAtoms.clear();
    Selection = NULL;
}

//------------------------------------------------------------------------------

bool CAmberMaskAtoms::SetMask(const CSmallString& mask)
{
    if( Topology == NULL ) {
        ES_ERROR("no topology is assigned with mask");
        return(false);
    }

    if( mask == NULL ) {
        ES_ERROR("mask is empty");
        return(false);
    }

    // update box information, required by distance operators
    if( Coordinates != NULL ) {
        Topology->BoxInfo.SetBoxDimmensions(Coordinates->GetBox());
        Topology->BoxInfo.SetBoxAngles(Coordinates->GetAngles());
        Topology->BoxInfo.UpdateBoxMatrices();
    }

    // remove previous selection -----------------------------
    Mask = mask;
    if( Selection != NULL ) delete Selection;
    Selection = NULL;
    SelectedAtoms.clear();

    // init mask parser
    init_mask();

    // parse mask
    if( parse_mask(Mask) != 0 ) {
        free_mask_tree();
        ES_ERROR("unable to parse mask");
        return(false);
    }

    // get top mask expression
    struct SExpression* p_top_expr = get_expression_tree();
    if( p_top_expr == NULL ) {
        ES_ERROR("top expression is NULL");
        return(false);
    }

    // init selection tree
    try {
        Selection = new CAmberMaskASelection(this);
    } catch(...) {
        free_mask_tree();
        ES_ERROR("unable to allocate root selection");
        return(false);
    }

    if( Selection->ExpandAndReduceTree(p_top_expr) == false ) {
        ES_ERROR("unable to expand and reduce expression tree");
        free_mask_tree();
        return(false);
    }

    // free parser data
    free_mask_tree();

    for(int i=0; i < GetNumberOfTopologyAtoms(); i++){
        CAmberAtom* p_atom = GetSelectedAtom(i);
        if( p_atom ){
            SelectedAtoms.push_back(p_atom);
        }
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberMaskAtoms::SetMaskFromFile(const CSmallString& maskfile)
{
    FILE* p_fin = fopen(maskfile,"r");
    if( p_fin == NULL ) {
        CSmallString error;
        error << "unable to open mask file '" << maskfile << "' (" << strerror(errno) << ")";
        ES_ERROR(error);
        return(false);
    }

    bool result = SetMaskFromFile(p_fin);

    fclose(p_fin);

    return(result);
}

//---------------------------------------------------------------------------

bool  CAmberMaskAtoms::SetMaskFromFile(FILE* p_file)
{
    CSmallString sbuffer;

    char buffer[512];
    bool end_loop = false;
    bool result = true;
    do {
        result = fgets(buffer,512,p_file) != NULL;
        if( result == true ) {
            int len = strlen(buffer);
            if( buffer[len-1] == '\n' ) {
                end_loop = true;
                buffer[len-1] = '\0';
            }
            sbuffer += buffer;
        }
    } while( (end_loop == false) && (result == true) );

    if( sbuffer == NULL ) {
        ES_ERROR("mask read from file is empty");
        return(false);
    }

    result = SetMask(sbuffer);

    return(result);
}

//---------------------------------------------------------------------------

const CSmallString& CAmberMaskAtoms::GetMask(void) const
{
    return(Mask);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberMaskAtoms::GetNumberOfTopologyAtoms(void)
{
    if( Topology == NULL ) return(0);
    return( Topology->AtomList.GetNumberOfAtoms() );
}

//------------------------------------------------------------------------------

int CAmberMaskAtoms::GetNumberOfSelectedAtoms(void)
{
    return( SelectedAtoms.size() );
}

//------------------------------------------------------------------------------

CAmberAtom* CAmberMaskAtoms::GetSelectedAtom(int index)
{
    if( Selection == NULL ) return(NULL);
    return( Selection->GetSelectedAtom(index) );
}

//------------------------------------------------------------------------------

bool CAmberMaskAtoms::IsAtomSelected(int index)
{
    if( Selection == NULL ) return(false);
    return( Selection->GetSelectedAtom(index) != NULL );
}

//------------------------------------------------------------------------------

CAmberAtom* CAmberMaskAtoms::GetSelectedAtomCondensed(int index)
{
    return( SelectedAtoms[index] );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberMaskAtoms::PrintInfo(bool short_info,FILE* fout)
{
    if( fout == NULL ) {
        fout = stdout;
    }

    fprintf(fout,"Atom based mask : %s\n",(const char*)Mask);
    fprintf(fout,"Number of atoms : %d\n",GetNumberOfSelectedAtoms());

    if( short_info == true ) {
        fprintf(fout,"\n");
        return;
    }

    CFortranIO  fortranio(fout);

    fortranio.SetFormat("3A25");
    fprintf(fout,"\n");
    fprintf(fout,"AtomID Name ResID Name   AtomID Name ResID Name   AtomID Name ResID Name\n");
    fprintf(fout,"======================   ======================   ======================\n");
    for(int i=0 ; i < GetNumberOfTopologyAtoms(); i++) {
        CAmberAtom* p_atm = GetSelectedAtom(i);
        if( p_atm != NULL ) {
            char buffer[50];
            sprintf(buffer,"%6d %4s %5d %4s   ",i+1,p_atm->GetName(),
                    p_atm->GetResidue()->GetIndex()+1,
                    p_atm->GetResidue()->GetName());
            fortranio.WriteString(buffer);
        }
    }
    fortranio.WriteEndOfSection();

    fprintf(fout,"\n");

    return;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


