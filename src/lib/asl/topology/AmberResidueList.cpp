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
//     This program is distributed","the hope that it will be useful,
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
#include <AmberResidueList.hpp>
#include <FortranIO.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberResidueList::CAmberResidueList(CAmberTopology* p_top)
{
    Topology = p_top;
    NRES = 0;
    NMXRS = 0;
    Residues = NULL;
    SortedResidues = NULL;
}

//------------------------------------------------------------------------------

CAmberResidueList::~CAmberResidueList(void)
{
    FreeFields();
    FreeSortedResidueList();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberResidueList::InitFields(int iNRES, int iNMXRS)
{
    FreeFields();
    FreeSortedResidueList();

    NRES = iNRES;
    NMXRS = iNMXRS;

    if( NRES > 0 ) Residues = new CAmberResidue[NRES];

    for(int i=0; i < NRES; i++) {
        Residues[i].Topology = Topology;
        Residues[i].Index = i;
    }
}

//------------------------------------------------------------------------------

void CAmberResidueList::FreeFields(void)
{
    if( Residues != NULL ) delete[] Residues;
    Residues = NULL;
    NRES = 0;
    NMXRS = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberResidueList::ReinitAtomResiduePointers(CAmberAtomList* p_atomlist)
{
    CAmberResidue* p_res = Residues;
    CAmberResidue* p_prev = NULL;

    for(int i=0; i<NRES; i++) {
        if( i != 0 ) {
            p_prev->NumOfAtoms = p_res->IPRES - p_prev->IPRES;
            CAmberAtom* p_atom = p_atomlist->GetAtom(p_prev->IPRES-1);
            for(int j=0; j<p_prev->NumOfAtoms; j++) {
                p_atom->Residue = p_prev;
                p_atom++;
            }
        }
        p_prev = p_res;
        p_res++;
    }

    p_prev->NumOfAtoms = p_atomlist->GetNumberOfAtoms() - p_prev->IPRES + 1;
    CAmberAtom* p_atom = p_atomlist->GetAtom(p_prev->IPRES-1);
    for(int j=0; j<p_prev->NumOfAtoms; j++) {
        p_atom->Residue = p_prev;
        p_atom++;
    }
}

//------------------------------------------------------------------------------

void CAmberResidueList::operator = (const CAmberResidueList& src)
{
//    void InitFields(int iNRES, int iNMXRS);
    InitFields(src.NRES,src.NMXRS);

    for(int i = 0; i < GetNumberOfResidues(); i++) {
        Residues[i] = src.Residues[i];
        Residues[i].Topology = Topology;
    }

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberResidueList::GetNumberOfResidues(void) const
{
    return(NRES);
}

//------------------------------------------------------------------------------

CAmberResidue* CAmberResidueList::GetResidue(int index) const
{
    return(&Residues[index]);
}

//---------------------------------------------------------------------------

bool CAmberResidueList::PrepareSortedResidues(void)
{
    FreeSortedResidueList(); // free previous list

    if( GetNumberOfResidues() > 0 ) SortedResidues = new CAmberResidue*[GetNumberOfResidues()];

    for(int i = 0; i < GetNumberOfResidues(); i++) {
        SortedResidues[i] = GetResidue(i);
    }

    qsort(SortedResidues,GetNumberOfResidues(),sizeof(CAmberResidue*),CompareItem);

    return(true);
}

//---------------------------------------------------------------------------

CAmberResidue* CAmberResidueList::GetSortedResidue(int index) const
{
    if( SortedResidues == NULL ) return(NULL);
    return(SortedResidues[index]);
}

//---------------------------------------------------------------------------

int CAmberResidueList::CompareItem(const void *p_left, const void *p_right)
{
    CAmberResidue* p_l = *((CAmberResidue**)p_left);
    CAmberResidue* p_r = *((CAmberResidue**)p_right);

    if( (p_l == NULL)||(p_r == NULL) ) return(0);
    return( strncmp(p_l->GetName(),p_r->GetName(),4) );
}

//------------------------------------------------------------------------------

bool CAmberResidueList::FreeSortedResidueList(void)
{
    if( SortedResidues != NULL ) delete[] SortedResidues;
    SortedResidues = NULL;
    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberResidueList::LoadResidueNames(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberResidue* p_res = Residues;

    for(int i=0; i<NRES; i++) {
        if( fortranio.ReadString(p_res->LABRES) == false ) {
            ES_ERROR("unable to load LABRES item");
            return(false);
        }
        p_res->Index = i;
        p_res++;
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberResidueList::LoadResidueIPRES(FILE* p_file,
        CAmberAtomList* p_atomlist,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberResidue* p_res = Residues;
    CAmberResidue* p_prev = NULL;

    for(int i=0; i<NRES; i++) {
        if( fortranio.ReadInt(p_res->IPRES) == false ) {
            ES_ERROR("unable to load IPRES item");
            return(false);
        }
        if( i != 0 ) {
            p_prev->NumOfAtoms = p_res->IPRES - p_prev->IPRES;
            CAmberAtom* p_atom = p_atomlist->GetAtom(p_prev->IPRES-1);
            for(int j=0; j<p_prev->NumOfAtoms; j++) {
                p_atom->Residue = p_prev;
                p_atom++;
            }
        }
        p_prev = p_res;
        p_res++;
    }

    p_prev->NumOfAtoms = p_atomlist->GetNumberOfAtoms() - p_prev->IPRES + 1;
    CAmberAtom* p_atom = p_atomlist->GetAtom(p_prev->IPRES-1);
    for(int j=0; j<p_prev->NumOfAtoms; j++) {
        p_atom->Residue = p_prev;
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberResidueList::LoadResiduePertNames(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberResidue* p_res = Residues;

    for(int i=0; i<NRES; i++) {
        if( fortranio.ReadString(p_res->PERRES) == false ) {
            ES_ERROR("unable to load PERRES item");
            return(false);
        }
        p_res++;
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberResidueList::SaveResidueNames(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberResidue* p_res = Residues;

    for(int i=0; i<NRES; i++) {
        if( fortranio.WriteString(p_res->LABRES) == false ) {
            ES_ERROR("unable to save LABRES item");
            return(false);
        }
        p_res++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberResidueList::SaveResidueIPRES(FILE* p_file,
        CAmberAtomList* p_atomlist,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberResidue* p_res = Residues;

    for(int i=0; i<NRES; i++) {
        if( fortranio.WriteInt(p_res->IPRES) == false ) {
            ES_ERROR("unable to save IPRES item");
            return(false);
        }
        p_res++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberResidueList::SaveResiduePertNames(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberResidue* p_res = Residues;

    for(int i=0; i<NRES; i++) {
        if( fortranio.WriteString(p_res->PERRES) == false ) {
            ES_ERROR("unable to save PERRES item");
            return(false);
        }
        p_res++;
    }
    fortranio.WriteEndOfSection();

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



