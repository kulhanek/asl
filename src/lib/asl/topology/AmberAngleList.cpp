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
#include <AmberAngle.hpp>
#include <AmberAngleList.hpp>
#include <FortranIO.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberAngleList::CAmberAngleList(void)
{
    NTHETH = 0;
    MTHETA = 0;
    NUMANG = 0;
    NTHETA = 0;
    NGPER = 0;
    MGPER = 0;

    AngleTypes = NULL;
    AngleWithoutHydrogens = NULL;
    AngleWithHydrogens = NULL;
    PerturbedAngles = NULL;
}

//------------------------------------------------------------------------------

CAmberAngleList::~CAmberAngleList(void)
{
    FreeFields();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberAngleList::InitFields(int iNTHETH, int iMTHETA,
        int iNUMANG, int iNTHETA, int iNGPER, int iMGPER)
{
    FreeFields();

    NTHETH = iNTHETH;
    MTHETA = iMTHETA;
    NUMANG = iNUMANG;
    NTHETA = iNTHETA;
    NGPER = iNGPER;
    MGPER = iMGPER;

    if( NUMANG > 0 ) AngleTypes = new CAmberAngleType[NUMANG];
    if( MTHETA > 0 ) AngleWithoutHydrogens = new CAmberAngle[MTHETA];
    if( NTHETH > 0 ) AngleWithHydrogens = new CAmberAngle[NTHETH];
    if( NGPER > 0 ) PerturbedAngles = new CAmberAngle[NGPER];
}

//------------------------------------------------------------------------------

void CAmberAngleList::FreeFields(void)
{
    if( AngleTypes != NULL ) delete[] AngleTypes;
    AngleTypes = NULL;

    if( AngleWithoutHydrogens != NULL ) delete[] AngleWithoutHydrogens;
    AngleWithoutHydrogens = NULL;

    if( AngleWithHydrogens != NULL ) delete[] AngleWithHydrogens;
    AngleWithHydrogens = NULL;

    if( PerturbedAngles != NULL ) delete[] PerturbedAngles;
    PerturbedAngles = NULL;

    NTHETH = 0;
    MTHETA = 0;
    NUMANG = 0;
    NTHETA = 0;
    NGPER = 0;
    MGPER = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberAngleList::GetNumberOfAngles(void) const
{
    return(NTHETH + MTHETA);
}

//------------------------------------------------------------------------------

int CAmberAngleList::GetNumberOfAnglesWithHydrogen(void) const
{
    return(NTHETH);
}

//------------------------------------------------------------------------------

int CAmberAngleList::GetNumberOfAnglesWithoutHydrogen(void) const
{
    return(MTHETA);
}

//------------------------------------------------------------------------------

int CAmberAngleList::GetNumberOfPerturbedAngles(void) const
{
    return(NGPER);
}

//------------------------------------------------------------------------------

int CAmberAngleList::GetNumberOfCompletelyPerturbedAngles(void) const
{
    return(MGPER);
}

//------------------------------------------------------------------------------

int CAmberAngleList::GetNumberOfAngleTypes(void) const
{
    return(NUMANG);
}

//------------------------------------------------------------------------------

int CAmberAngleList::GetNTHETA(void) const
{
    return(NTHETA);
}

//------------------------------------------------------------------------------

CAmberAngle* CAmberAngleList::GetAngleWithHydrogen(int index) const
{
    return(&AngleWithHydrogens[index]);
}

//------------------------------------------------------------------------------

CAmberAngle* CAmberAngleList::GetAngleWithoutHydrogen(int index) const
{
    return(&AngleWithoutHydrogens[index]);
}

//------------------------------------------------------------------------------

CAmberAngle* CAmberAngleList::GetAngle(int index) const
{
    if( index < MTHETA ) {
        return(&AngleWithoutHydrogens[index]);
    } else {
        return(&AngleWithHydrogens[index-MTHETA]);
    }
}

//------------------------------------------------------------------------------

CAmberAngle* CAmberAngleList::GetPerturbedAngle(int index) const
{
    return(&PerturbedAngles[index]);
}

//------------------------------------------------------------------------------

CAmberAngleType* CAmberAngleList::GetAngleType(int index) const
{
    return(&AngleTypes[index]);
}

//------------------------------------------------------------------------------

void CAmberAngleList::operator = (const CAmberAngleList& src)
{
//    void InitFields(int iNTHETH, int iMTHETA, int iNUMANG,
//                                 int iNTHETA, int iNGPER, int iMGPER);
    InitFields(src.NTHETH,src.MTHETA,src.NUMANG,src.NTHETA,src.NGPER,src.MGPER);

//    int NTHETH;   //NTHETH : number of angles containing hydrogen
//    int MTHETA;   //MTHETA : number of angles not containing hydrogen
//    int NUMANG;   //NUMANG : number of unique angle types
//    int NTHETA;   //NTHETA : MTHETA + number of constraint angles
//    int NGPER;    //NGPER  : number of angles to be perturbed
//    int MGPER;    //MGPER  : number of angles with atoms completely in perturbed group

    for(int i=0; i<NUMANG; i++){
        AngleTypes[i] = src.AngleTypes[i];
    }
    for(int i=0; i<NTHETH; i++){
        AngleWithHydrogens[i] = src.AngleWithHydrogens[i];
    }
    for(int i=0; i<MTHETA; i++){
        AngleWithoutHydrogens[i] = src.AngleWithoutHydrogens[i];
    }
    for(int i=0; i<NGPER; i++){
        PerturbedAngles[i] = src.PerturbedAngles[i];
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberAngleList::LoadAngleTK(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngleType* p_type = AngleTypes;

    for(int i=0; i<NUMANG; i++) {
        if( fortranio.ReadReal(p_type->TK) == false ) {
            ES_ERROR("unable to load TK item");
            return(false);
        }
        p_type++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAngleList::LoadAngleTEQ(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngleType* p_type = AngleTypes;

    for(int i=0; i<NUMANG; i++) {
        if( fortranio.ReadReal(p_type->TEQ) == false ) {
            ES_ERROR("unable to load TEQ item");
            return(false);
        }
        p_type++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAngleList::LoadAnglesWithHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngle* p_angle = AngleWithHydrogens;

    for(int i=0; i<NTHETH; i++) {
        if( fortranio.ReadInt(p_angle->IT) == false ) {
            ES_ERROR("unable to load IT item");
            return(false);
        }
        if( fortranio.ReadInt(p_angle->JT) == false ) {
            ES_ERROR("unable to load JT item");
            return(false);
        }
        if( fortranio.ReadInt(p_angle->KT) == false ) {
            ES_ERROR("unable to load KT item");
            return(false);
        }
        if( fortranio.ReadInt(p_angle->ICT) == false ) {
            ES_ERROR("unable to load ICT item");
            return(false);
        }

        // reindex
        p_angle->IT = p_angle->IT/3;
        p_angle->JT = p_angle->JT/3;
        p_angle->KT = p_angle->KT/3;
        p_angle->ICT = p_angle->ICT - 1;

        p_angle++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAngleList::LoadAnglesWithoutHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngle* p_angle = AngleWithoutHydrogens;

    for(int i=0; i<MTHETA; i++) {
        if( fortranio.ReadInt(p_angle->IT) == false ) {
            ES_ERROR("unable to load IT item");
            return(false);
        }
        if( fortranio.ReadInt(p_angle->JT) == false ) {
            ES_ERROR("unable to load JT item");
            return(false);
        }
        if( fortranio.ReadInt(p_angle->KT) == false ) {
            ES_ERROR("unable to load KT item");
            return(false);
        }
        if( fortranio.ReadInt(p_angle->ICT) == false ) {
            ES_ERROR("unable to load ICT item");
            return(false);
        }

        // reindex
        p_angle->IT = p_angle->IT/3;
        p_angle->JT = p_angle->JT/3;
        p_angle->KT = p_angle->KT/3;
        p_angle->ICT = p_angle->ICT - 1;

        p_angle++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAngleList::LoadPerturbedAngles(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngle* p_angle = PerturbedAngles;

    for(int i=0; i<NGPER; i++) {
        if( fortranio.ReadInt(p_angle->IT) == false ) {
            ES_ERROR("unable to load IT item");
            return(false);
        }
        if( fortranio.ReadInt(p_angle->JT) == false ) {
            ES_ERROR("unable to load JT item");
            return(false);
        }
        if( fortranio.ReadInt(p_angle->KT) == false ) {
            ES_ERROR("unable to load KT item");
            return(false);
        }

        // reindex
        p_angle->IT = p_angle->IT/3;
        p_angle->JT = p_angle->JT/3;
        p_angle->KT = p_angle->KT/3;

        p_angle++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAngleList::LoadPerturbedAngleTypeIndexes(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngle* p_angle = PerturbedAngles;

    for(int i=0; i<NGPER; i++) {
        if( fortranio.ReadInt(p_angle->ICT) == false ) {
            ES_ERROR("unable to load ICT item");
            return(false);
        }
        // reindex
        p_angle->ICT = p_angle->ICT - 1;

        p_angle++;
    }

    p_angle = PerturbedAngles;

    for(int i=0; i<NGPER; i++) {
        if( fortranio.ReadInt(p_angle->PCT) == false ) {
            ES_ERROR("unable to load PCT item");
            return(false);
        }
        // reindex
        p_angle->PCT = p_angle->PCT - 1;

        p_angle++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAngleList::SaveAngleTK(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngleType* p_type = AngleTypes;

    for(int i=0; i<NUMANG; i++) {
        if( fortranio.WriteReal(p_type->TK) == false ) {
            ES_ERROR("unable to save TK item");
            return(false);
        }
        p_type++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAngleList::SaveAngleTEQ(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngleType* p_type = AngleTypes;

    for(int i=0; i<NUMANG; i++) {
        if( fortranio.WriteReal(p_type->TEQ) == false ) {
            ES_ERROR("unable to save TEQ item");
            return(false);
        }
        p_type++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAngleList::SaveAnglesWithHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngle* p_angle = AngleWithHydrogens;

    for(int i=0; i<NTHETH; i++) {
        if( fortranio.WriteInt(p_angle->IT*3) == false ) {
            ES_ERROR("unable to save IT item");
            return(false);
        }
        if( fortranio.WriteInt(p_angle->JT*3) == false ) {
            ES_ERROR("unable to save JT item");
            return(false);
        }
        if( fortranio.WriteInt(p_angle->KT*3) == false ) {
            ES_ERROR("unable to save KT item");
            return(false);
        }
        if( fortranio.WriteInt(p_angle->ICT+1) == false ) {
            ES_ERROR("unable to save ICT item");
            return(false);
        }

        p_angle++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAngleList::SaveAnglesWithoutHydrogens(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngle* p_angle = AngleWithoutHydrogens;

    for(int i=0; i<MTHETA; i++) {
        if( fortranio.WriteInt(p_angle->IT*3) == false ) {
            ES_ERROR("unable to save IT item");
            return(false);
        }
        if( fortranio.WriteInt(p_angle->JT*3) == false ) {
            ES_ERROR("unable to save JT item");
            return(false);
        }
        if( fortranio.WriteInt(p_angle->KT*3) == false ) {
            ES_ERROR("unable to save KT item");
            return(false);
        }
        if( fortranio.WriteInt(p_angle->ICT+1) == false ) {
            ES_ERROR("unable to save ICT item");
            return(false);
        }

        p_angle++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAngleList::SavePerturbedAngles(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngle* p_angle = PerturbedAngles;

    for(int i=0; i<NGPER; i++) {
        if( fortranio.WriteInt(p_angle->IT*3) == false ) {
            ES_ERROR("unable to save IT item");
            return(false);
        }
        if( fortranio.WriteInt(p_angle->JT*3) == false ) {
            ES_ERROR("unable to save JT item");
            return(false);
        }
        if( fortranio.WriteInt(p_angle->KT*3) == false ) {
            ES_ERROR("unable to save KT item");
            return(false);
        }

        p_angle++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAngleList::SavePerturbedAngleTypeIndexes(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAngle* p_angle = PerturbedAngles;

    for(int i=0; i<NGPER; i++) {
        if( fortranio.WriteInt(p_angle->ICT + 1) == false ) {
            ES_ERROR("unable to save ICT item");
            return(false);
        }
        p_angle++;
    }

    p_angle = PerturbedAngles;
    for(int i=0; i<NGPER; i++) {
        if( fortranio.WriteInt(p_angle->PCT + 1) == false ) {
            ES_ERROR("unable to save PCT item");
            return(false);
        }
        p_angle++;
    }

    fortranio.WriteEndOfSection();
    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


