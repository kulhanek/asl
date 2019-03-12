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
#include <AmberDihedralList.hpp>
#include <FortranIO.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberDihedralList::CAmberDihedralList(void)
{
    NPHIH = 0;
    MPHIA = 0;
    NPTRA = 0;
    NPHIA = 0;
    NDPER = 0;
    MDPER = 0;

    DihedralTypes = NULL;
    DihedralWithoutHydrogens = NULL;
    DihedralWithHydrogens = NULL;
    PerturbedDihedrals = NULL;

    SCEEFactorsLoaded = false;
    SCNBFactorsLoaded = false;
}

//------------------------------------------------------------------------------

CAmberDihedralList::~CAmberDihedralList(void)
{
    FreeFields();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberDihedralList::InitFields(int iNPHIH, int iMPHIA,int iNPTRA,
                                    int iNPHIA, int iNDPER, int iMDPER)
{
    FreeFields();

    NPHIH = iNPHIH;
    MPHIA = iMPHIA;
    NPTRA = iNPTRA;
    NPHIA = iNPHIA;
    NDPER = iNDPER;
    MDPER = iMDPER;

    if( NPTRA > 0 ) DihedralTypes = new CAmberDihedralType[NPTRA];
    if( MPHIA > 0 ) DihedralWithoutHydrogens = new CAmberDihedral[MPHIA];
    if( NPHIH > 0 ) DihedralWithHydrogens = new CAmberDihedral[NPHIH];
    if( NDPER > 0 ) PerturbedDihedrals = new CAmberDihedral[NDPER];
}

//------------------------------------------------------------------------------

void CAmberDihedralList::FreeFields(void)
{
    if( DihedralTypes != NULL ) delete[] DihedralTypes;
    DihedralTypes = NULL;

    if( DihedralWithoutHydrogens != NULL ) delete[] DihedralWithoutHydrogens;
    DihedralWithoutHydrogens = NULL;

    if( DihedralWithHydrogens != NULL ) delete[] DihedralWithHydrogens;
    DihedralWithHydrogens = NULL;

    if( PerturbedDihedrals != NULL ) delete[] PerturbedDihedrals;
    PerturbedDihedrals = NULL;

    NPHIH = 0;
    MPHIA = 0;
    NPTRA = 0;
    NPHIA = 0;
    NDPER = 0;
    MDPER = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberDihedralList::GetNumberOfDihedrals(void) const
{
    return(NPHIH + MPHIA);
}

//------------------------------------------------------------------------------

int CAmberDihedralList::GetNumberOfDihedralsWithHydrogen(void) const
{
    return(NPHIH);
}

//------------------------------------------------------------------------------

int CAmberDihedralList::GetNumberOfDihedralsWithoutHydrogen(void) const
{
    return(MPHIA);
}

//------------------------------------------------------------------------------

int CAmberDihedralList::GetNumberOfPerturbedDihedrals(void) const
{
    return(NDPER);
}

//------------------------------------------------------------------------------

int CAmberDihedralList::GetNumberOfCompletelyPerturbedDihedrals(void) const
{
    return(MDPER);
}

//------------------------------------------------------------------------------

int CAmberDihedralList::GetNumberOfDihedralTypes(void) const
{
    return(NPTRA);
}

//------------------------------------------------------------------------------

int CAmberDihedralList::GetNPHIA(void) const
{
    return(NPHIA);
}

//------------------------------------------------------------------------------

CAmberDihedral* CAmberDihedralList::GetDihedralWithHydrogen(int index) const
{
    return(&DihedralWithHydrogens[index]);
}

//------------------------------------------------------------------------------

CAmberDihedral* CAmberDihedralList::GetDihedralWithoutHydrogen(int index) const
{
    return(&DihedralWithoutHydrogens[index]);
}

//------------------------------------------------------------------------------

CAmberDihedral* CAmberDihedralList::GetDihedral(int index) const
{
    if( index < MPHIA ) {
        return(&DihedralWithoutHydrogens[index]);

    } else {
        return(&DihedralWithHydrogens[index-MPHIA]);
    }
}

//------------------------------------------------------------------------------

CAmberDihedral* CAmberDihedralList::GetPerturbedDihedral(int index) const
{
    return(&PerturbedDihedrals[index]);
}

//------------------------------------------------------------------------------

CAmberDihedralType* CAmberDihedralList::GetDihedralType(int index) const
{
    return(&DihedralTypes[index]);
}

//------------------------------------------------------------------------------

void CAmberDihedralList::operator = (const CAmberDihedralList& src)
{
//    void InitFields(int iNPHIH, int iMPHIA, int iNPTRA,
//                    int iNPHIA, int iNDPER, int iMDPER);
    InitFields(src.NPHIH,src.MPHIA,src.NPTRA,src.NPHIA,src.NDPER,src.MDPER);

//    int NPHIH;    //NPHIH  : number of dihedrals containing hydrogen
//    int MPHIA;    //MPHIA  : number of dihedrals not containing hydrogen
//    int NPTRA;    //NPTRA  : number of unique dihedral types
//    int NPHIA;    //NPHIA  : MPHIA + number of constraint dihedrals
//    int NDPER;    //NDPER  : number of dihedrals to be perturbed
//    int MDPER;    //MDPER  : number of dihedrals with atoms completely in perturbed groups

    for(int i=0; i<NPTRA; i++) {
        DihedralTypes[i] = src.DihedralTypes[i];
    }

    for(int i=0; i<MPHIA; i++) {
        DihedralWithoutHydrogens[i] = src.DihedralWithoutHydrogens[i];
    }

    for(int i=0; i<NPHIH; i++) {
        DihedralWithHydrogens[i] = src.DihedralWithHydrogens[i];
    }

    for(int i=0; i<NDPER; i++) {
        PerturbedDihedrals[i] = src.PerturbedDihedrals[i];
    }
}

//------------------------------------------------------------------------------

void CAmberDihedralList::RemoveIllegalDihedrals(void)
{
    CAmberDihedral* OldDihedralWithHydrogens = DihedralWithHydrogens;
    int             OldNPHIH = NPHIH;
    CAmberDihedral* OldDihedralWithoutHydrogens = DihedralWithoutHydrogens;
    int             OldMPHIA = MPHIA;

    // recalculate number of angles
    NPHIH = 0;
    for(int i=0; i < OldNPHIH; i++) {
        if( OldDihedralWithHydrogens[i].GetICP() >= 0 ) NPHIH++;
    }
    MPHIA = 0;
    for(int i=0; i < OldMPHIA; i++) {
        if( OldDihedralWithoutHydrogens[i].GetICP() >= 0 ) MPHIA++;
    }

    DihedralWithHydrogens = NULL;
    if( NPHIH > 0 ) DihedralWithHydrogens = new CAmberDihedral[NPHIH];
    DihedralWithoutHydrogens = NULL;
    if( MPHIA > 0 ) DihedralWithoutHydrogens = new CAmberDihedral[MPHIA];

    int j = 0;
    for(int i=0; i < OldNPHIH; i++) {
        if( OldDihedralWithHydrogens[i].GetICP() >= 0 ){
            DihedralWithHydrogens[j] = OldDihedralWithHydrogens[i];
            j++;
        }
    }

    j = 0;
    for(int i=0; i < OldMPHIA; i++) {
        if( OldDihedralWithoutHydrogens[i].GetICP() >= 0 ){
            DihedralWithoutHydrogens[j] = OldDihedralWithoutHydrogens[i];
            j++;
        }
    }

    // clean old data
    if( OldDihedralWithHydrogens )     delete[] OldDihedralWithHydrogens;
    if( OldDihedralWithoutHydrogens )  delete[] OldDihedralWithoutHydrogens;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberDihedralList::LoadDihedralPK(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedralType* p_type = DihedralTypes;

    for(int i=0; i<NPTRA; i++) {
        if( fortranio.ReadReal(p_type->PK) == false ) {
            ES_ERROR("unable to load TEQ item");
            return(false);
        }
        p_type++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::LoadDihedralPN(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedralType* p_type = DihedralTypes;

    for(int i=0; i<NPTRA; i++) {
        if( fortranio.ReadReal(p_type->PN) == false ) {
            ES_ERROR("unable to load PN item");
            return(false);
        }
        p_type++;
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::LoadDihedralPHASE(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedralType* p_type = DihedralTypes;

    for(int i=0; i<NPTRA; i++) {
        if( fortranio.ReadReal(p_type->PHASE) == false ) {
            ES_ERROR("unable to load PHASE item");
            return(false);
        }
        p_type++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::LoadDihedralSCEE(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedralType* p_type = DihedralTypes;

    for(int i=0; i<NPTRA; i++) {
        if( fortranio.ReadReal(p_type->SCEE_SCALE) == false ) {
            ES_ERROR("unable to load SCEE_SCALE item");
            return(false);
        }
        p_type++;
    }
    SCEEFactorsLoaded = true;
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::LoadDihedralSCNB(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedralType* p_type = DihedralTypes;

    for(int i=0; i<NPTRA; i++) {
        if( fortranio.ReadReal(p_type->SCNB_SCALE) == false ) {
            ES_ERROR("unable to load SCNB_SCALE item");
            return(false);
        }
        p_type++;
    }
    SCNBFactorsLoaded = true;
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::LoadDihedralsWithHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedral* p_dihedral = DihedralWithHydrogens;

    for(int i=0; i<NPHIH; i++) {
        if( fortranio.ReadInt(p_dihedral->IP) == false ) {
            ES_ERROR("unable to load IP item");
            return(false);
        }
        if( fortranio.ReadInt(p_dihedral->JP) == false ) {
            ES_ERROR("unable to load JP item");
            return(false);
        }
        if( fortranio.ReadInt(p_dihedral->KP) == false ) {
            ES_ERROR("unable to load KP item");
            return(false);
        }
        if( fortranio.ReadInt(p_dihedral->LP) == false ) {
            ES_ERROR("unable to load LP item");
            return(false);
        }
        if( fortranio.ReadInt(p_dihedral->ICP) == false ) {
            ES_ERROR("unable to load ICP item");
            return(false);
        }

        // reindex
        p_dihedral->IP = p_dihedral->IP/3;
        p_dihedral->JP = p_dihedral->JP/3;
        p_dihedral->KP = p_dihedral->KP/3;
        p_dihedral->LP = p_dihedral->LP/3;

        if( (p_dihedral->KP < 0) && (p_dihedral->LP < 0) ) {
            p_dihedral->Type = -2;
            p_dihedral->LP *= -1;
            p_dihedral->KP *= -1;
        }
        if( p_dihedral->KP < 0 ) {
            p_dihedral->Type = -1;
            p_dihedral->KP *= -1;
        }
        if( p_dihedral->LP < 0 ) {
            p_dihedral->Type = 1;
            p_dihedral->LP *= -1;
        }
        p_dihedral->ICP = p_dihedral->ICP - 1;

        p_dihedral++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::LoadDihedralsWithoutHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedral* p_dihedral = DihedralWithoutHydrogens;

    for(int i=0; i<MPHIA; i++) {
        if( fortranio.ReadInt(p_dihedral->IP) == false ) {
            ES_ERROR("unable to load IP item");
            return(false);
        }
        if( fortranio.ReadInt(p_dihedral->JP) == false ) {
            ES_ERROR("unable to load JP item");
            return(false);
        }
        if( fortranio.ReadInt(p_dihedral->KP) == false ) {
            ES_ERROR("unable to load KP item");
            return(false);
        }
        if( fortranio.ReadInt(p_dihedral->LP) == false ) {
            ES_ERROR("unable to load LP item");
            return(false);
        }
        if( fortranio.ReadInt(p_dihedral->ICP) == false ) {
            ES_ERROR("unable to load ICP item");
            return(false);
        }

        // reindex
        p_dihedral->IP = p_dihedral->IP/3;
        p_dihedral->JP = p_dihedral->JP/3;
        p_dihedral->KP = p_dihedral->KP/3;
        p_dihedral->LP = p_dihedral->LP/3;
        if( (p_dihedral->KP < 0) && (p_dihedral->LP < 0) ) {
            p_dihedral->Type = -2;
            p_dihedral->LP *= -1;
            p_dihedral->KP *= -1;
        }
        if( p_dihedral->KP < 0 ) {
            p_dihedral->Type = -1;
            p_dihedral->KP *= -1;
        }
        if( p_dihedral->LP < 0 ) {
            p_dihedral->Type = 1;
            p_dihedral->LP *= -1;
        }
        p_dihedral->ICP = p_dihedral->ICP - 1;

        p_dihedral++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::LoadPerturbedDihedrals(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedral* p_dihedral = PerturbedDihedrals;

    for(int i=0; i<NDPER; i++) {
        if( fortranio.ReadInt(p_dihedral->IP) == false ) {
            ES_ERROR("unable to load IP item");
            return(false);
        }
        if( fortranio.ReadInt(p_dihedral->JP) == false ) {
            ES_ERROR("unable to load JP item");
            return(false);
        }
        if( fortranio.ReadInt(p_dihedral->KP) == false ) {
            ES_ERROR("unable to load KP item");
            return(false);
        }
        if( fortranio.ReadInt(p_dihedral->LP) == false ) {
            ES_ERROR("unable to load LP item");
            return(false);
        }

        // reindex
        p_dihedral->IP = p_dihedral->IP/3;
        p_dihedral->JP = p_dihedral->JP/3;
        p_dihedral->KP = p_dihedral->KP/3;
        p_dihedral->LP = p_dihedral->LP/3;
        if( (p_dihedral->KP < 0) && (p_dihedral->LP < 0) ) {
            p_dihedral->Type = -2;
            p_dihedral->LP *= -1;
            p_dihedral->KP *= -1;
        }

        if( p_dihedral->KP < 0 ) {
            p_dihedral->Type = -1;
            p_dihedral->KP *= -1;
        }
        if( p_dihedral->LP < 0 ) {
            p_dihedral->Type = 1;
            p_dihedral->LP *= -1;
        }

        p_dihedral++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::LoadPerturbedDihedralTypeIndexes(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedral* p_dihedral = PerturbedDihedrals;

    for(int i=0; i<NDPER; i++) {
        if( fortranio.ReadInt(p_dihedral->ICP) == false ) {
            ES_ERROR("unable to load ICP item");
            return(false);
        }

        // reindex
        p_dihedral->ICP = p_dihedral->ICP  - 1;

        p_dihedral++;
    }

    p_dihedral = PerturbedDihedrals;
    for(int i=0; i<NDPER; i++) {
        if( fortranio.ReadInt(p_dihedral->PCP) == false ) {
            ES_ERROR("unable to load PCP item");
            return(false);
        }

        // reindex
        p_dihedral->PCP = p_dihedral->PCP  - 1;

        p_dihedral++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::SaveDihedralPK(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedralType* p_type = DihedralTypes;

    for(int i=0; i<NPTRA; i++) {
        if( fortranio.WriteReal(p_type->PK) == false ) {
            ES_ERROR("unable to save TEQ item");
            return(false);
        }
        p_type++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::SaveDihedralPN(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedralType* p_type = DihedralTypes;

    for(int i=0; i<NPTRA; i++) {
        if( fortranio.WriteReal(p_type->PN) == false ) {
            ES_ERROR("unable to save PN item");
            return(false);
        }
        p_type++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::SaveDihedralPHASE(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedralType* p_type = DihedralTypes;

    for(int i=0; i<NPTRA; i++) {
        if( fortranio.WriteReal(p_type->PHASE) == false ) {
            ES_ERROR("unable to save PHASE item");
            return(false);
        }
        p_type++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::SaveDihedralSCEE(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedralType* p_type = DihedralTypes;

    for(int i=0; i<NPTRA; i++) {
        if( fortranio.WriteReal(p_type->SCEE_SCALE) == false ) {
            ES_ERROR("unable to save SCEE_SCALE item");
            return(false);
        }
        p_type++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::SaveDihedralSCNB(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedralType* p_type = DihedralTypes;

    for(int i=0; i<NPTRA; i++) {
        if( fortranio.WriteReal(p_type->SCNB_SCALE) == false ) {
            ES_ERROR("unable to save SCNB_SCALE item");
            return(false);
        }
        p_type++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::SaveDihedralsWithHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedral* p_dihedral = DihedralWithHydrogens;
    int value;

    for(int i=0; i<NPHIH; i++) {
        if( fortranio.WriteInt(p_dihedral->IP*3) == false ) {
            ES_ERROR("unable to save IP item");
            return(false);
        }
        if( fortranio.WriteInt(p_dihedral->JP*3) == false ) {
            ES_ERROR("unable to save JP item");
            return(false);
        }

        value = p_dihedral->KP;
        if( (p_dihedral->Type == -1) || (p_dihedral->Type == -2) ) value *= -1;
        if( fortranio.WriteInt(value*3) == false ) {
            ES_ERROR("unable to save KP item");
            return(false);
        }

        value = p_dihedral->LP;
        if( (p_dihedral->Type == 1) || (p_dihedral->Type == -2) ) value *= -1;
        if( fortranio.WriteInt(value*3) == false ) {
            ES_ERROR("unable to save LP item");
            return(false);
        }

        if( fortranio.WriteInt(p_dihedral->ICP+1) == false ) {
            ES_ERROR("unable to save ICP item");
            return(false);
        }
        p_dihedral++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberDihedralList::SaveDihedralsWithoutHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedral* p_dihedral = DihedralWithoutHydrogens;

    int value;

    for(int i=0; i<MPHIA; i++) {
        if( fortranio.WriteInt(p_dihedral->IP*3) == false ) {
            ES_ERROR("unable to save IP item");
            return(false);
        }
        if( fortranio.WriteInt(p_dihedral->JP*3) == false ) {
            ES_ERROR("unable to save JP item");
            return(false);
        }

        value = p_dihedral->KP;
        if( (p_dihedral->Type == -1) || (p_dihedral->Type == -2) ) value *= -1;
        if( fortranio.WriteInt(value*3) == false ) {
            ES_ERROR("unable to save KP item");
            return(false);
        }

        value = p_dihedral->LP;
        if( (p_dihedral->Type == 1) || (p_dihedral->Type == -2) ) value *= -1;
        if( fortranio.WriteInt(value*3) == false ) {
            ES_ERROR("unable to save LP item");
            return(false);
        }
        if( fortranio.WriteInt(p_dihedral->ICP+1) == false ) {
            ES_ERROR("unable to save ICP item");
            return(false);
        }

        p_dihedral++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberDihedralList::SavePerturbedDihedrals(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedral* p_dihedral = PerturbedDihedrals;

    int value;

    for(int i=0; i<NDPER; i++) {
        if( fortranio.WriteInt(p_dihedral->IP*3) == false ) {
            ES_ERROR("unable to save IP item");
            return(false);
        }
        if( fortranio.WriteInt(p_dihedral->JP*3) == false ) {
            ES_ERROR("unable to save JP item");
            return(false);
        }

        value = p_dihedral->KP;
        if( p_dihedral->Type == -1 ) value *= -1;
        if( fortranio.WriteInt(p_dihedral->KP) == false ) {
            ES_ERROR("unable to save KP item");
            return(false);
        }

        value = p_dihedral->LP;
        if( p_dihedral->Type == 1 ) value *= -1;
        if( fortranio.WriteInt(p_dihedral->LP) == false ) {
            ES_ERROR("unable to save LP item");
            return(false);
        }

        p_dihedral++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//---------------------------------------------------------------------------

bool CAmberDihedralList::SavePerturbedDihedralTypeIndexes(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberDihedral* p_dihedral = PerturbedDihedrals;

    for(int i=0; i<NDPER; i++) {
        if( fortranio.WriteInt(p_dihedral->ICP + 1) == false ) {
            ES_ERROR("unable to save ICP item");
            return(false);
        }
        p_dihedral++;
    }

    p_dihedral = PerturbedDihedrals;
    for(int i=0; i<NDPER; i++) {
        if( fortranio.WriteInt(p_dihedral->PCP + 1) == false ) {
            ES_ERROR("unable to save ICP item");
            return(false);
        }
        p_dihedral++;
    }

    fortranio.WriteEndOfSection();

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

