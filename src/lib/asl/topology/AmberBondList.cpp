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
#include <AmberBondList.hpp>
#include <FortranIO.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberBondList::CAmberBondList(void)
{
    NBONH = 0;
    MBONA = 0;
    NUMBND = 0;
    NBONA = 0;
    NBPER = 0;
    MBPER = 0;

    BondTypes = NULL;
    BondsWithoutHydrogens = NULL;
    BondsWithHydrogens = NULL;
    PerturbedBonds = NULL;
}

//------------------------------------------------------------------------------

CAmberBondList::~CAmberBondList(void)
{
    FreeFields();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberBondList::InitFields(int iNBONH,int iMBONA,int iNUMBND,
        int iNBONA,int iNBPER,int iMBPER)
{
    FreeFields();

    NBONH = iNBONH;
    MBONA = iMBONA;
    NUMBND = iNUMBND;
    NBONA = iNBONA;
    NBPER = iNBPER;
    MBPER = iMBPER;

    if( NUMBND > 0 ) BondTypes = new CAmberBondType[NUMBND];
    if( MBONA > 0 ) BondsWithoutHydrogens = new CAmberBond[MBONA];
    if( NBONH > 0 ) BondsWithHydrogens = new CAmberBond[NBONH];
    if( NBPER > 0 ) PerturbedBonds = new CAmberBond[NBPER];
}

//------------------------------------------------------------------------------

void CAmberBondList::FreeFields(void)
{
    if( BondTypes != NULL ) delete[] BondTypes;
    BondTypes = NULL;

    if( BondsWithoutHydrogens != NULL ) delete[] BondsWithoutHydrogens;
    BondsWithoutHydrogens = NULL;

    if( BondsWithHydrogens != NULL ) delete[] BondsWithHydrogens;
    BondsWithHydrogens = NULL;

    if( PerturbedBonds != NULL ) delete[] PerturbedBonds;
    PerturbedBonds = NULL;

    NBONH = 0;
    MBONA = 0;
    NUMBND = 0;
    NBONA = 0;
    NBPER = 0;
    MBPER = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberBondList::GetNumberOfBonds(void) const
{
    return(MBONA + NBONH);
}

//------------------------------------------------------------------------------

int CAmberBondList::GetNumberOfBondsWithHydrogen(void) const
{
    return(NBONH);
}

//------------------------------------------------------------------------------

int CAmberBondList::GetNumberOfBondsWithoutHydrogen(void) const
{
    return(MBONA);
}

//------------------------------------------------------------------------------

int CAmberBondList::GetNumberOfPerturbedBonds(void) const
{
    return(NBPER);
}

//------------------------------------------------------------------------------

int CAmberBondList::GetNumberOfCompletelyPerturbedBonds(void) const
{
    return(MBPER);
}

//------------------------------------------------------------------------------

int CAmberBondList::GetNumberOfBondTypes(void) const
{
    return(NUMBND);
}

//------------------------------------------------------------------------------

int CAmberBondList::GetNBONA(void) const
{
    return(NBONA);
}

//------------------------------------------------------------------------------

CAmberBond* CAmberBondList::GetBondWithHydrogen(int index) const
{
    return(&BondsWithHydrogens[index]);
}

//------------------------------------------------------------------------------

CAmberBond* CAmberBondList::GetBondWithoutHydrogen(int index) const
{
    return(&BondsWithoutHydrogens[index]);
}

//------------------------------------------------------------------------------

CAmberBond* CAmberBondList::GetBond(int index) const
{
    if( index < MBONA ) {
        return(&BondsWithoutHydrogens[index]);
    } else {
        return(&BondsWithHydrogens[index-MBONA]);
    }
}

//------------------------------------------------------------------------------

CAmberBond* CAmberBondList::GetPerturbedBond(int index) const
{
    return(&PerturbedBonds[index]);
}

//------------------------------------------------------------------------------

CAmberBondType* CAmberBondList::GetBondType(int index) const
{
    return(&BondTypes[index]);
}

//------------------------------------------------------------------------------

void CAmberBondList::operator = (const CAmberBondList& src)
{
//    void InitFields(int iNBONH,int iMBONA,int iNUMBND,int iNBONA,
//                    int iNBPER,int iMBPER);
    InitFields(src.NBONH,src.MBONA,src.NUMBND,src.NBONA,src.NBPER,src.MBPER);

//    int     NBONH;    //NBONH  : number of bonds containing hydrogen
//    int     MBONA;    //MBONA  : number of bonds not containing hydrogen
//    int     NUMBND;   //NUMBND : number of unique bond types
//    int     NBONA;    //NBONA  : MBONA + number of constraint bonds
//    int     NBPER;    //NBPER  : number of bonds to be perturbed
//    int     MBPER;    //MBPER  : number of bonds with atoms completely in perturbed group

    for(int i=0; i<NUMBND; i++) {
        BondTypes[i] = src.BondTypes[i];
    }
    for(int i=0; i<MBONA; i++) {
        BondsWithoutHydrogens[i] = src.BondsWithoutHydrogens[i];
    }
    for(int i=0; i<NBONH; i++) {
        BondsWithHydrogens[i] = src.BondsWithHydrogens[i];
    }
    for(int i=0; i<NBPER; i++) {
        PerturbedBonds[i] = src.PerturbedBonds[i];
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberBondList::LoadBondRK(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBondType* p_type = BondTypes;

    for(int i=0; i<NUMBND; i++) {
        if( fortranio.ReadReal(p_type->RK) == false ) {
            ES_ERROR("unable to load RK item");
            return(false);
        }
        p_type++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBondList::LoadBondREQ(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBondType* p_type = BondTypes;

    for(int i=0; i<NUMBND; i++) {
        if( fortranio.ReadReal(p_type->REQ) == false ) {
            ES_ERROR("unable to load REQ item");
            return(false);
        }
        p_type++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBondList::LoadBondsWithHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBond* p_bond = BondsWithHydrogens;

    for(int i=0; i<NBONH; i++) {
        if( fortranio.ReadInt(p_bond->IB) == false ) {
            ES_ERROR("unable to load IB item");
            return(false);
        }
        if( fortranio.ReadInt(p_bond->JB) == false ) {
            ES_ERROR("unable to load JB item");
            return(false);
        }
        if( fortranio.ReadInt(p_bond->ICB) == false ) {
            ES_ERROR("unable to load ICB item");
            return(false);
        }

        // reindex
        p_bond->IB = p_bond->IB/3;
        p_bond->JB = p_bond->JB/3;
        p_bond->ICB = p_bond->ICB - 1;

        p_bond++;
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBondList::LoadBondsWithoutHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBond* p_bond = BondsWithoutHydrogens;

    for(int i=0; i<MBONA; i++) {
        if( fortranio.ReadInt(p_bond->IB) == false ) {
            ES_ERROR("unable to load IB item");
            return(false);
        }
        if( fortranio.ReadInt(p_bond->JB) == false ) {
            ES_ERROR("unable to load JB item");
            return(false);
        }
        if( fortranio.ReadInt(p_bond->ICB) == false ) {
            ES_ERROR("unable to load ICB item");
            return(false);
        }

        // reindex
        p_bond->IB = p_bond->IB/3;
        p_bond->JB = p_bond->JB/3;
        p_bond->ICB = p_bond->ICB - 1;

        p_bond++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBondList::LoadPerturbedBonds(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBond* p_bond = PerturbedBonds;

    for(int i=0; i<NBPER; i++) {
        if( fortranio.ReadInt(p_bond->IB) == false ) {
            ES_ERROR("unable to load IB item");
            return(false);
        }
        if( fortranio.ReadInt(p_bond->JB) == false ) {
            ES_ERROR("unable to load JB item");
            return(false);
        }
        // reindex
        p_bond->IB = p_bond->IB/3;
        p_bond->JB = p_bond->JB/3;

        p_bond++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBondList::LoadPerturbedBondTypeIndexes(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBond* p_bond = PerturbedBonds;

    for(int i=0; i<NBPER; i++) {
        if( fortranio.ReadInt(p_bond->ICB) == false ) {
            ES_ERROR("unable to load ICB item");
            return(false);
        }
        // reindex
        p_bond->ICB = p_bond->ICB - 1;

        p_bond++;
    }

    p_bond = PerturbedBonds;

    for(int i=0; i<NBPER; i++) {
        if( fortranio.ReadInt(p_bond->PCB) == false ) {
            ES_ERROR("unable to load PCB item");
            return(false);
        }
        // reindex
        p_bond->PCB = p_bond->PCB - 1;

        p_bond++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBondList::SaveBondRK(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBondType* p_type = BondTypes;

    for(int i=0; i<NUMBND; i++) {
        if( fortranio.WriteReal(p_type->RK) == false ) {
            ES_ERROR("unable to save RK item");
            return(false);
        }
        p_type++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBondList::SaveBondREQ(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBondType* p_type = BondTypes;

    for(int i=0; i<NUMBND; i++) {
        if( fortranio.WriteReal(p_type->REQ) == false ) {
            ES_ERROR("unable to save REQ item");
            return(false);
        }
        p_type++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBondList::SaveBondsWithHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBond* p_bond = BondsWithHydrogens;

    for(int i=0; i<NBONH; i++) {
        if( fortranio.WriteInt(p_bond->IB*3) == false ) {
            ES_ERROR("unable to save IB item");
            return(false);
        }
        if( fortranio.WriteInt(p_bond->JB*3) == false ) {
            ES_ERROR("unable to saveJB item");
            return(false);
        }
        if( fortranio.WriteInt(p_bond->ICB+1) == false ) {
            ES_ERROR("unable to saveICB item");
            return(false);
        }

        p_bond++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBondList::SaveBondsWithoutHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBond* p_bond = BondsWithoutHydrogens;

    for(int i=0; i<MBONA; i++) {
        if( fortranio.WriteInt(p_bond->IB*3) == false ) {
            ES_ERROR("unable to save IB item");
            return(false);
        }
        if( fortranio.WriteInt(p_bond->JB*3) == false ) {
            ES_ERROR("unable to save JB item");
            return(false);
        }
        if( fortranio.WriteInt(p_bond->ICB+1) == false ) {
            ES_ERROR("unable to save ICB item");
            return(false);
        }

        p_bond++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberBondList::SavePerturbedBonds(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBond* p_bond = PerturbedBonds;

    for(int i=0; i<NBPER; i++) {
        if( fortranio.WriteInt(p_bond->IB*3) == false ) {
            ES_ERROR("unable to save IB item");
            return(false);
        }
        if( fortranio.WriteInt(p_bond->JB*3) == false ) {
            ES_ERROR("unable to save JB item");
            return(false);
        }
        p_bond++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberBondList::SavePerturbedBondTypeIndexes(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberBond* p_bond = PerturbedBonds;

    for(int i=0; i<NBPER; i++) {
        if( fortranio.WriteInt(p_bond->ICB+1) == false ) {
            ES_ERROR("unable to save ICB item");
            return(false);
        }
        p_bond++;
    }

    p_bond = PerturbedBonds;
    for(int i=0; i<NBPER; i++) {
        if( fortranio.WriteInt(p_bond->PCB+1) == false ) {
            ES_ERROR("unable to save PCB item");
            return(false);
        }
        p_bond++;
    }

    fortranio.WriteEndOfSection();
    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

