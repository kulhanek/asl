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
#include <AmberNonBondedList.hpp>
#include <FortranIO.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberNonBondedList::CAmberNonBondedList(void)
{
    ICO = NULL;
    CN1 = NULL;
    CN2 = NULL;
    NATEX = NULL;
    ASOL = NULL;
    BSOL = NULL;
    HBCUT = NULL;
    SOLTY = NULL;
    NTYPES = 0;
    NEXT = 0;
    NATYP = 0;
    NPHB = 0;
}

//------------------------------------------------------------------------------

CAmberNonBondedList::~CAmberNonBondedList(void)
{
    FreeFields();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberNonBondedList::InitFields(int iNTYPES, int iNEXT,
        int iNATYP, int iNPHB)
{
    FreeFields(); // destroy previous data

    NTYPES = iNTYPES;
    NEXT = iNEXT;
    NATYP = iNATYP;
    NPHB = iNPHB;

    // allocate neccessary fields
    if( NTYPES*NTYPES > 0 ) ICO = new int[NTYPES*NTYPES];

    if( NTYPES*(NTYPES+1)/2 > 0 ) CN1 = new double[NTYPES*(NTYPES+1)/2];
    for(int i=0; i < NTYPES*(NTYPES+1)/2; i++ ) CN1[i] = 0.0;

    if( NTYPES*(NTYPES+1)/2 > 0 ) CN2 = new double[NTYPES*(NTYPES+1)/2];
    for(int i=0; i < NTYPES*(NTYPES+1)/2; i++ ) CN2[i] = 0.0;

    if( NEXT > 0 ) NATEX = new int[NEXT];

    if( NPHB > 0 ) ASOL = new double[NPHB];
    for(int i=0; i < NPHB; i++ ) ASOL[i] = 0.0;

    if( NPHB > 0 ) BSOL = new double[NPHB];
    for(int i=0; i < NPHB; i++ ) BSOL[i] = 0.0;

    if( NPHB > 0 ) HBCUT = new double[NPHB];
    for(int i=0; i < NPHB; i++ ) HBCUT[i] = 0.0;

    if( NATYP > 0 ) SOLTY = new double[NATYP];
    for(int i=0; i < NATYP; i++ ) SOLTY[i] = 0.0;
}

//------------------------------------------------------------------------------

void CAmberNonBondedList::FreeFields(void)
{
    if( ICO != NULL ) delete[] ICO;
    if( CN1 != NULL ) delete[] CN1;
    if( CN2 != NULL ) delete[] CN2;
    if( NATEX != NULL ) delete[] NATEX;
    if( ASOL != NULL ) delete[] ASOL;
    if( BSOL != NULL ) delete[] BSOL;
    if( HBCUT != NULL ) delete[] HBCUT;
    if( SOLTY != NULL ) delete[] SOLTY;
    ICO = NULL;
    CN1 = NULL;
    CN2 = NULL;
    NATEX = NULL;
    ASOL = NULL;
    BSOL = NULL;
    HBCUT = NULL;
    SOLTY = NULL;
    NTYPES = 0;
    NEXT = 0;
    NATYP = 0;
    NPHB = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberNonBondedList::GetNumberOfTypes(void)
{
    return(NTYPES);
}

//------------------------------------------------------------------------------

int CAmberNonBondedList::GetNumberOfExcludedAtoms(void)
{
    return(NEXT);
}

//------------------------------------------------------------------------------

int CAmberNonBondedList::GetNATEX(int index)
{
    return(NATEX[index]);
}

//------------------------------------------------------------------------------

void CAmberNonBondedList::SetNATEX(int index,int value)
{
    NATEX[index] = value;
}

//------------------------------------------------------------------------------

int  CAmberNonBondedList::GetNonBondedType(int index)
{
    int ico = ICO[index];
    if( ico != 0 ) ico = ico / abs(ico); // sign(ico)
    return(ico);
}

//------------------------------------------------------------------------------

int  CAmberNonBondedList::GetNonBondedType(CAmberAtom* p_at1,CAmberAtom* p_at2,bool pert)
{
    if( (p_at1 == NULL) || (p_at2 == NULL) ) return(0);

    int ico = ICO[NTYPES*(p_at1->GetIAC(pert)-1)+p_at2->GetIAC(pert)-1];
    if( ico != 0 ) ico = ico / abs(ico);  // sign(ico)

    return(ico);
}

//------------------------------------------------------------------------------

double CAmberNonBondedList::GetAParam(CAmberAtom* p_at1,CAmberAtom* p_at2,bool pert)
{
    if( (p_at1 == NULL) || (p_at2 == NULL) ) return(0);

    int ico = ICO[NTYPES*(p_at1->GetIAC(pert)-1)+p_at2->GetIAC(pert)-1];
    int type = ico;
    if( ico != 0 ) ico = abs(ico) - 1;
    if( type > 0 ) return(CN1[ico]);
    else   return(ASOL[ico]);
}

//------------------------------------------------------------------------------

double CAmberNonBondedList::GetBParam(CAmberAtom* p_at1,CAmberAtom* p_at2,bool pert)
{
    if( (p_at1 == NULL) || (p_at2 == NULL) ) return(0);

    int ico = ICO[NTYPES*(p_at1->GetIAC(pert)-1)+p_at2->GetIAC(pert)-1];
    int type = ico;
    if( ico != 0 ) ico = abs(ico) - 1;
    if( type > 0 ) return(CN2[ico]);
    else   return(BSOL[ico]);
}

//------------------------------------------------------------------------------

int CAmberNonBondedList::GetICOIndex(int index)
{
    return(ICO[index]);
}

//------------------------------------------------------------------------------

void CAmberNonBondedList::SetICOIndex(int index,int value)
{
    ICO[index] = value;
}

//------------------------------------------------------------------------------

int CAmberNonBondedList::GetICOIndex(CAmberAtom* p_at1,CAmberAtom* p_at2,bool pert)
{
    if( (p_at1 == NULL) || (p_at2 == NULL) ) return(0);

    int ico = ICO[NTYPES*(p_at1->GetIAC(pert)-1)+p_at2->GetIAC(pert)-1];

    return(ico);
}

//------------------------------------------------------------------------------

double CAmberNonBondedList::GetAParam(int icoindex)
{
    if( icoindex == 0 ) return(0);
    int type = icoindex;
    icoindex = abs(icoindex) - 1;
    if( type > 0 ) return(CN1[icoindex]);
    else   return(ASOL[icoindex]);
}

//------------------------------------------------------------------------------

double CAmberNonBondedList::GetBParam(int icoindex)
{
    if( icoindex == 0 ) return(0);
    int type = icoindex;
    icoindex = abs(icoindex) - 1;
    if( type > 0 ) return(CN2[icoindex]);
    else   return(BSOL[icoindex]);
}

//------------------------------------------------------------------------------

void CAmberNonBondedList::SetAParam(double value,int icoindex)
{
    if( icoindex == 0 ) return;
    int type = icoindex;
    icoindex = abs(icoindex) - 1;
    if( type > 0 ) CN1[icoindex] = value;
    else   ASOL[icoindex] = value;
}

//------------------------------------------------------------------------------

void CAmberNonBondedList::SetBParam(double value,int icoindex)
{
    if( icoindex == 0 ) return;
    int type = icoindex;
    icoindex = abs(icoindex) - 1;
    if( type > 0 ) CN2[icoindex] = value;
    else   BSOL[icoindex] = value;
}

//------------------------------------------------------------------------------

void CAmberNonBondedList::operator = (const CAmberNonBondedList& src)
{
//void InitFields(int iNTYPES, int iNEXT, int iNATYP, int iNPHB);
    InitFields(src.NTYPES,src.NEXT,src.NATYP,src.NPHB);

    for(int i=0; i<NTYPES*NTYPES; i++) {
        ICO[i] = src.ICO[i];
    }

    for(int i=0; i<NATYP; i++) {
        SOLTY[i] = src.SOLTY[i];
    }

    for(int i=0; i<NTYPES*(NTYPES+1)/2; i++) {
        CN1[i] = src.CN1[i];
        CN1[i] = src.CN1[i];
    }

    for(int i=0; i<NEXT; i++) {
        NATEX[i] = src.NATEX[i];
    }

    for(int i=0; i<NPHB; i++) {
        ASOL[i] = src.ASOL[i];
        BSOL[i] = src.BSOL[i];
        HBCUT[i] = src.HBCUT[i];
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberNonBondedList::LoadICOs(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    int* p_item = ICO;

    for(int i=0; i<NTYPES*NTYPES; i++) {
        if( fortranio.ReadInt(*p_item) == false ) {
            ES_ERROR("unable to load ICO item");
            return(false);
        }
        p_item++;
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::LoadSOLTY(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = SOLTY;

    for(int i=0; i<NATYP; i++) {
        if( fortranio.ReadReal(*p_item) == false ) {
            ES_ERROR("unable to load SOLTY item");
            return(false);
        }
        p_item++;
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::LoadCN1(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = CN1;

    for(int i=0; i<NTYPES*(NTYPES+1)/2; i++) {
        if( fortranio.ReadReal(*p_item) == false ) {
            ES_ERROR("unable to load CN1 item");
            return(false);
        }
        p_item++;
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::LoadCN2(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = CN2;

    for(int i=0; i<NTYPES*(NTYPES+1)/2; i++) {
        if( fortranio.ReadReal(*p_item) == false ) {
            ES_ERROR("unable to load CN2 item");
            return(false);
        }
        p_item++;
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::LoadNATEX(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    int* p_item = NATEX;

    for(int i=0; i<NEXT; i++) {
        if( fortranio.ReadInt(*p_item) == false ) {
            ES_ERROR("unable to load NATEX item");
            return(false);
        }
        *p_item = *p_item - 1;
        p_item++;
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::LoadASOL(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = ASOL;

    for(int i=0; i<NPHB; i++) {
        if( fortranio.ReadReal(*p_item) == false ) {
            ES_ERROR("unable to load ASOL item");
            return(false);
        }
        p_item++;
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::LoadBSOL(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = BSOL;

    for(int i=0; i<NPHB; i++) {
        if( fortranio.ReadReal(*p_item) == false ) {
            ES_ERROR("unable to load BSOL item");
            return(false);
        }
        p_item++;
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::LoadHBCUT(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = HBCUT;

    for(int i=0; i<NPHB; i++) {
        if( fortranio.ReadReal(*p_item) == false ) {
            ES_ERROR("unable to load HBCUT item");
            return(false);
        }
        p_item++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberNonBondedList::SaveICOs(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    int* p_item = ICO;

    for(int i=0; i<NTYPES*NTYPES; i++) {
        if( fortranio.WriteInt(*p_item) == false ) {
            ES_ERROR("unable to save item");
            return(false);
        }
        p_item++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::SaveSOLTY(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = SOLTY;

    for(int i=0; i<NATYP; i++) {
        if( fortranio.WriteReal(*p_item) == false ) {
            ES_ERROR("unable to save SOLTY item");
            return(false);
        }
        p_item++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::SaveCN1(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = CN1;

    for(int i=0; i<NTYPES*(NTYPES+1)/2; i++) {
        if( fortranio.WriteReal(*p_item) == false ) {
            ES_ERROR("unable to save CN1 item");
            return(false);
        }
        p_item++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::SaveCN2(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = CN2;

    for(int i=0; i<NTYPES*(NTYPES+1)/2; i++) {
        if( fortranio.WriteReal(*p_item) == false ) {
            ES_ERROR("unable to save CN2 item");
            return(false);
        }
        p_item++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::SaveNATEX(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    int* p_item = NATEX;

    for(int i=0; i<NEXT; i++) {
        if( fortranio.WriteInt((*p_item) + 1) == false ) {
            ES_ERROR("unable to save NATEX item");
            return(false);
        }
        p_item++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::SaveASOL(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = ASOL;

    for(int i=0; i<NPHB; i++) {
        if( fortranio.WriteReal(*p_item) == false ) {
            ES_ERROR("unable to save ASOL item");
            return(false);
        }
        p_item++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::SaveBSOL(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = BSOL;

    for(int i=0; i<NPHB; i++) {
        if( fortranio.WriteReal(*p_item) == false ) {
            ES_ERROR("unable to save BSOL item");
            return(false);
        }
        p_item++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberNonBondedList::SaveHBCUT(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    double* p_item = HBCUT;

    for(int i=0; i<NPHB; i++) {
        if( fortranio.WriteReal(*p_item) == false ) {
            ES_ERROR("unable to save HBCUT item");
            return(false);
        }
        p_item++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------
//##############################################################################
//------------------------------------------------------------------------------



