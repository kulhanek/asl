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
#include <AmberAtomList.hpp>
#include <FortranIO.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberAtomList::CAmberAtomList(void)
{
    NATOM = 0;
    IFPERT = 0;
    IFPOL = 0;
    Atoms =  NULL;
    AtomicNumberLoaded = false;
}

//------------------------------------------------------------------------------

CAmberAtomList::~CAmberAtomList(void)
{
    FreeFields();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberAtomList::InitFields(int iNATOM,int iIFPERT,int iIFPOL)
{
    FreeFields();

    NATOM = iNATOM;
    IFPERT = iIFPERT;
    IFPOL = iIFPOL;

    if( NATOM > 0 ) Atoms = new CAmberAtom[NATOM];

    for(int i=0; i < NATOM; i++) Atoms[i].AtomIndex = i;
}

//------------------------------------------------------------------------------

void CAmberAtomList::FreeFields(void)
{
    if( Atoms != NULL ) delete[] Atoms;
    Atoms = NULL;

    NATOM = 0;
    IFPERT = 0;
    IFPOL = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberAtomList::GetNumberOfAtoms(void) const
{
    return(NATOM);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::HasPertInfo(void) const
{
    return(IFPERT == 1);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::HasPolInfo(void) const
{
    return(IFPOL == 1);
}

//------------------------------------------------------------------------------

CAmberAtom* CAmberAtomList::GetAtom(int index) const
{
    return(&Atoms[index]);
}

//------------------------------------------------------------------------------

const CSmallString& CAmberAtomList::GetRadiusSet(void) const
{
    return(RadiusSet);
}

//------------------------------------------------------------------------------

void CAmberAtomList::SetRadiusSet(const CSmallString& set_name)
{
    RadiusSet = set_name;
}

//------------------------------------------------------------------------------

void CAmberAtomList::operator = (const CAmberAtomList& src)
{
//    void InitFields(int iNATOM,int iIFPERT,int iIFPOL);
    InitFields(src.NATOM,src.IFPERT,src.IFPOL);

    // copy name of radius set
    RadiusSet = src.RadiusSet;

    // do we have atomic numbers?
    AtomicNumberLoaded = src.AtomicNumberLoaded;

    // copy atoms
    for(int i=0; i<NATOM; i++){
        Atoms[i] = src.Atoms[i];
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberAtomList::LoadAtomNames(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadString(p_atom->IGRAPH) == false ) {
            ES_ERROR("unable to load IGRAPH item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomCharges(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadReal(p_atom->CHRG) == false ) {
            ES_ERROR("unable to load CHRG item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomAtomicNumbers(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadInt(p_atom->ATOMIC_NUMBER) == false ) {
            ES_ERROR("unable to load ATOMIC_NUMBER item");
            return(false);
        }
        p_atom++;
    }
    AtomicNumberLoaded = true;
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomMasses(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadReal(p_atom->AMASS) == false ) {
            ES_ERROR("unable to load AMASS item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomIACs(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadInt(p_atom->IAC) == false ) {
            ES_ERROR("unable to load IAC item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomNUMEXs(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadInt(p_atom->NUMEX) == false ) {
            ES_ERROR("unable to load NUMEX item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomIPol(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);
    return(fortranio.ReadInt(IFPOL));
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomPol(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadReal(p_atom->ATPOL) == false ) {
            ES_ERROR("unable to load ATPOL item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomNames(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteString(p_atom->IGRAPH) == false ) {
            ES_ERROR("unable to save IGRAPH item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomCharges(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteReal(p_atom->CHRG) == false ) {
            ES_ERROR("unable to save CHRG item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomAtomicNumbers(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteInt(p_atom->ATOMIC_NUMBER) == false ) {
            ES_ERROR("unable to save ATOMIC_NUMBER item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomMasses(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteReal(p_atom->AMASS) == false ) {
            ES_ERROR("unable to save AMASS item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomIACs(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteInt(p_atom->IAC) == false ) {
            ES_ERROR("unable to save AMASS item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomNUMEXs(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteInt(p_atom->NUMEX) == false ) {
            ES_ERROR("unable to save NUMEX item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomPol(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteReal(p_atom->ATPOL) == false ) {
            ES_ERROR("unable to save ATPOL item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomISYMBL(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadString(p_atom->ISYMBL) == false ) {
            ES_ERROR("unable to load ISYMBL item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomITREE(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadString(p_atom->ITREE) == false ) {
            ES_ERROR("unable to load ISYMBL item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomJOIN(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadInt(p_atom->JOIN) == false ) {
            ES_ERROR("unable to load JOIN item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomIROTAT(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadInt(p_atom->IROTAT) == false ) {
            ES_ERROR("unable to load IROTAT item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomISYMBL(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteString(p_atom->ISYMBL) == false ) {
            ES_ERROR("unable to save ISYMBL item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomITREE(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteString(p_atom->ITREE) == false ) {
            ES_ERROR("unable to save ISYMBL item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomJOIN(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteInt(p_atom->JOIN) == false ) {
            ES_ERROR("unable to save JOIN item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomIROTAT(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteInt(p_atom->IROTAT) == false ) {
            ES_ERROR("unable to save IROTAT item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadPertAtomNames(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadString(p_atom->IGRPER) == false ) {
            ES_ERROR("unable to load IGRPER item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadPertAtomISYMBL(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadString(p_atom->ISMPER) == false ) {
            ES_ERROR("unable to load ISMPER item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadPertAtomCharges(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadReal(p_atom->CGPER) == false ) {
            ES_ERROR("unable to load CGPER item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadPertAtomPertFlag(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadInt(p_atom->IAPER) == false ) {
            ES_ERROR("unable to load IAPER item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadPertAtomIAC(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadInt(p_atom->IACPER) == false ) {
            ES_ERROR("unable to load IACPER item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadPertAtomPol(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadReal(p_atom->ATPOL1) == false ) {
            ES_ERROR("unable to load ATPOL1 item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadPertAtomALMPER(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadReal(p_atom->ALMPER) == false ) {
            ES_ERROR("unable to load ALMPER item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SavePertAtomNames(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteString(p_atom->IGRPER) == false ) {
            ES_ERROR("unable to save IGRPER item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SavePertAtomISYMBL(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteString(p_atom->ISMPER) == false ) {
            ES_ERROR("unable to save ISMPER item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SavePertAtomCharges(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteReal(p_atom->CGPER) == false ) {
            ES_ERROR("unable to save CGPER item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SavePertAtomPertFlag(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteInt(p_atom->IAPER) == false ) {
            ES_ERROR("unable to save IAPER item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SavePertAtomIAC(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteInt(p_atom->IACPER) == false ) {
            ES_ERROR("unable to save IACPER item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SavePertAtomPol(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteReal(p_atom->ATPOL1) == false ) {
            ES_ERROR("unable to save ATPOL1 item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SavePertAtomALMPER(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteReal(p_atom->ALMPER) == false ) {
            ES_ERROR("unable to save ALMPER item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomRadiusSet(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);
    return(fortranio.ReadString(RadiusSet));
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomRadii(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadReal(p_atom->RADIUS) == false ) {
            ES_ERROR("unable to load RADIUS item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::LoadAtomScreen(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.ReadReal(p_atom->SCREEN) == false ) {
            ES_ERROR("unable to load SCREEN item");
            return(false);
        }
        p_atom++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomRadiusSet(FILE* p_file,const char* p_format)
{
    if( RadiusSet == NULL ) return(true);

    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    fortranio.WriteString(RadiusSet);

    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomRadii(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteReal(p_atom->RADIUS) == false ) {
            ES_ERROR("unable to save RADIUS item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberAtomList::SaveAtomScreen(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    CAmberAtom* p_atom = Atoms;

    for(int i=0; i<NATOM; i++) {
        if( fortranio.WriteReal(p_atom->SCREEN) == false ) {
            ES_ERROR("unable to save SCREEN item");
            return(false);
        }
        p_atom++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




