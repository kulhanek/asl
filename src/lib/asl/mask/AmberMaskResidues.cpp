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

#include <AmberMaskResidues.hpp>
#include <AmberTopology.hpp>
#include <AmberRestart.hpp>
#include <AmberMaskRSelection.hpp>
#include <FortranIO.hpp>
#include <ErrorSystem.hpp>

#include "maskparser/AmberMaskParser.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberMaskResidues::CAmberMaskResidues(void)
{
    Topology = NULL;
    Coordinates = NULL;
    Selection = NULL;
}

//------------------------------------------------------------------------------

CAmberMaskResidues::~CAmberMaskResidues(void)
{
    if( Selection != NULL ) delete Selection;
    Selection = NULL;
    Topology = NULL;
    Coordinates = NULL;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberMaskResidues::AssignTopology(CAmberTopology* p_top)
{
    if( Selection != NULL ) delete Selection;
    SelectedResidues.clear();
    Topology = p_top;
    Selection = NULL;
    Coordinates = NULL;
    Mask = NULL;
    return(true);
}

//------------------------------------------------------------------------------

CAmberTopology* CAmberMaskResidues::GetTopology(void) const
{
    return(Topology);
}

//------------------------------------------------------------------------------

bool CAmberMaskResidues::AssignCoordinates(CAmberRestart* p_crd)
{
    Coordinates = p_crd;
    return(true);
}

//------------------------------------------------------------------------------

CAmberRestart* CAmberMaskResidues::GetCoordinates(void) const
{
    return(Coordinates);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberMaskResidues::SelectAllResidues(void)
{
    return(SetMask(":*"));
}

//------------------------------------------------------------------------------

void CAmberMaskResidues::Reset(void)
{
    if( Selection != NULL ) delete Selection;
    SelectedResidues.clear();
    Selection = NULL;
}

//------------------------------------------------------------------------------

bool CAmberMaskResidues::SetMask(const CSmallString& mask)
{
    if( Topology == NULL ) {
        ES_ERROR("No topology is assigned with mask!");
        return(false);
    }

    if( mask == NULL ) {
        ES_ERROR("Mask is empty!");
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
    SelectedResidues.clear();

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
        Selection = new CAmberMaskRSelection(this);
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

    for(int i=0; i < GetNumberOfTopologyResidues(); i++){
        CAmberResidue* p_res = GetSelectedResidue(i);
        if( p_res ){
            SelectedResidues.push_back(p_res);
        }
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberMaskResidues::SetMaskFromFile(const CSmallString& maskfile)
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

bool  CAmberMaskResidues::SetMaskFromFile(FILE* p_file)
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

const CSmallString& CAmberMaskResidues::GetMask(void) const
{
    return(Mask);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberMaskResidues::GetNumberOfTopologyResidues(void)
{
    if( Topology == NULL ) return(0);
    return( Topology->ResidueList.GetNumberOfResidues() );
}

//------------------------------------------------------------------------------

int CAmberMaskResidues::GetNumberOfSelectedResidues(void)
{
    return( SelectedResidues.size() );
}

//------------------------------------------------------------------------------

CAmberResidue* CAmberMaskResidues::GetSelectedResidue(int index)
{
    if( Selection == NULL ) return(NULL);
    return( Selection->GetSelectedResidue(index) );
}

//------------------------------------------------------------------------------

bool CAmberMaskResidues::IsResidueSelected(int index)
{
    if( Selection == NULL ) return(false);
    return( Selection->GetSelectedResidue(index) != NULL );
}

//------------------------------------------------------------------------------

CAmberResidue* CAmberMaskResidues::GetSelectedResidueCondensed(int index)
{
    return( SelectedResidues[index] );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberMaskResidues::PrintInfo(bool short_info,FILE* fout)
{
    if( fout == NULL ) {
        fout = stdout;
    }

    fprintf(fout,"Residue based mask : %s\n",(const char*)Mask);
    fprintf(fout,"Number of residues : %d\n",GetNumberOfSelectedResidues());

    if( short_info == true ) {
        fprintf(fout,"\n");
        return;
    }

    CFortranIO  fortranio(fout);

    fortranio.SetFormat("3A12");
    fprintf(fout,"\n");
    fprintf(fout,"ResID Name   ResID Name   ResID Name\n");
    fprintf(fout,"==========   ==========   ==========\n");
    for(int i=0 ; i < GetNumberOfTopologyResidues(); i++) {
        CAmberResidue* p_res = GetSelectedResidue(i);
        if( p_res != NULL ) {
            char buffer[50];
            sprintf(buffer,"%5d %4s   ",
                    p_res->GetIndex()+1,
                    p_res->GetName());
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


