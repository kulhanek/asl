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

#include <AmberRestart.hpp>
#include <AmberTopology.hpp>
#include <FortranIO.hpp>
#include <string.h>
#include <ErrorSystem.hpp>
#include <errno.h>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>

//---------------------------------------------------------------------------

CPoint CAmberRestart::zero;  // returned value when fields are not allocated

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberRestart::CAmberRestart(void)
{
    Topology = NULL;
    Positions = NULL;
    Velocities = NULL;
    VelocitiesLoaded = false;
    Title=NULL;
    Time=0;
    NumberOfAtoms=0;
}

//---------------------------------------------------------------------------

CAmberRestart::~CAmberRestart(void)
{
    Release();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberRestart::AssignTopology(CAmberTopology* p_topology)
{
    Topology = p_topology;
    Release();
}

//------------------------------------------------------------------------------

int CAmberRestart::GetNumberOfAtoms(void) const
{
    if( Positions ) return(NumberOfAtoms);
    return(0);
}

//------------------------------------------------------------------------------

CAmberTopology* CAmberRestart::GetTopology(void)
{
    return(Topology);
}

//---------------------------------------------------------------------------

bool CAmberRestart::Create(void)
{
    Release();
    if( Topology == NULL ) return(false);
    if( Topology->AtomList.GetNumberOfAtoms() <= 0 ) return(false);

    NumberOfAtoms = Topology->AtomList.GetNumberOfAtoms();

    Positions = new CPoint[NumberOfAtoms];
    if( Positions == NULL ) return(false);
    memset(Positions,0,NumberOfAtoms*sizeof(CPoint));

    Velocities = new CPoint[NumberOfAtoms];
    if( Velocities == NULL ) return(false);
    memset(Velocities,0,NumberOfAtoms*sizeof(CPoint));

    Box.x = 0.0;
    Box.y = 0.0;
    Box.z = 0.0;

    Box1.x = 0.0;
    Box1.y = 0.0;
    Box1.z = 0.0;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberRestart::Release(void)
{
    if( Positions != NULL ) delete[]  Positions;
    Positions = NULL;
    if( Velocities != NULL ) delete[] Velocities;
    Velocities = NULL;
    Box.x = 0.0;
    Box.y = 0.0;
    Box.z = 0.0;
    Box1.x = 0.0;
    Box1.y = 0.0;
    Box1.z = 0.0;
    VelocitiesLoaded = false;
    Title=NULL;
    Time=0;
    NumberOfAtoms=0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberRestart::Load(const CSmallString& name,bool allow_stdin)
{
    FILE* fin;

    if( (allow_stdin == true) && (name == "-") ) {
        fin = stdin;
    } else {
        fin = fopen(name,"r");
        if( fin == NULL ) {
            CSmallString error;
            error << "unable to open restart file '" << name << "' ("
                  << strerror(errno) << ")";
            ES_ERROR(error);
            return(false);
        }
    }

    bool  result = Load(fin);

    if( ! ((allow_stdin == true) && (name == "-")) ) {
        fclose(fin);
    }

    return(result);
}

//---------------------------------------------------------------------------

bool CAmberRestart::Save(const CSmallString& name,bool allow_stdout)
{
    FILE* fout;

    if( (allow_stdout == true) && (name == "-") ) {
        fout = stdout;
    } else {
        fout = fopen(name,"w");
        if( fout == NULL ) {
            CSmallString error;
            error << "unable to open restart file '" << name << "' ("
                  << strerror(errno) << ")";
            ES_ERROR(error);
            return(false);
        }
    }

    bool  result = Save(fout);

    if( ! ((allow_stdout == true) && (name == "-")) ) {
        fclose(fout);
    }

    return(result);
}

//---------------------------------------------------------------------------

bool CAmberRestart::Load(FILE *fin)
{
    if( fin == NULL ) return(false);

    Release();

    if( Topology == NULL ) {
        ES_ERROR("topology is not specified");
        return(false);
    }
    if( Topology->AtomList.GetNumberOfAtoms() == 0 ) {
        ES_ERROR("topology does not have any atom - is it loaded?");
        return(false);
    }

    if( Create() == false ) {
        ES_ERROR("unable to allocate memory for restart file");
        return(false);
    }

    // load title
    char buffer[100];  // be sure that \n is also loaded
    if( fgets(buffer,100,fin) == NULL ) {
        ES_ERROR("unable to load restart title");
        return(false);
    }

    buffer[80] = '\0';
    // remove trailing /n character
    buffer[strlen(buffer)-1] = '\0';

    Title = buffer;

    CFortranIO fortranio(fin,true);

    // this is old format for AMBER6
    fortranio.SetFormat("1I5");

    if( fortranio.ReadInt(NumberOfAtoms) == false ) {
        // it was not passed
        // then try 1I6
        fortranio.ReturnBack();
        fortranio.ChangeFormat("1I6");
        if( fortranio.ReadInt(NumberOfAtoms) == false ) {
            ES_ERROR("unable to load number of atoms from restart file");
            return(false);
        }
    }

    if( Topology->AtomList.GetNumberOfAtoms() != NumberOfAtoms ) {
        // it was not passed
        // then try 1I6
        fortranio.ReturnBack();
        fortranio.ChangeFormat("1I6");
        if( fortranio.ReadInt(NumberOfAtoms) == false ) {
            ES_ERROR("unable to load number of atoms from restart file");
            return(false);
        }
        if( Topology->AtomList.GetNumberOfAtoms() != NumberOfAtoms ) {
            CSmallString error;
            error << "number of atoms in topology " << Topology->AtomList.GetNumberOfAtoms() <<
                  " does not match the number of atoms in restart file " << NumberOfAtoms;
            ES_ERROR(error);
            return(false);
        }
    }

    fortranio.ChangeFormat("1E15.7");
    if( fortranio.ReadReal(Time) == false ) {
        Time = 0.0;
        // ES_ERROR("Unable to load time position from restart file !\n");
        // return(false);
        // restart file from xleap does not contain time
    }

    fortranio.SetFormat("6F12.7");
    bool result = true;
    for(int i=0; i < NumberOfAtoms; i++) {
        result &= fortranio.ReadReal(Positions[i].x);
        result &= fortranio.ReadReal(Positions[i].y);
        result &= fortranio.ReadReal(Positions[i].z);
        if( result == false ) {
            CSmallString error;
            error << "unable to load coordinates of atom " << i;
            ES_ERROR(error);
            return(false);
        }
    }
    fortranio.SetFormat("6F12.7");
    VelocitiesLoaded = true;
    int loaded_v = 0;
    for(int i=0; i < NumberOfAtoms; i++) {
        result &= fortranio.ReadReal(Velocities[i].x);
        if( result == true ) loaded_v++;
        result &= fortranio.ReadReal(Velocities[i].y);
        if( result == true ) loaded_v++;
        result &= fortranio.ReadReal(Velocities[i].z);
        if( result == true ) loaded_v++;
        if( result == false ) {
            VelocitiesLoaded = false;
            break;
        }
    }
    if( (VelocitiesLoaded == true) && (Topology->BoxInfo.GetType() != AMBER_BOX_NONE) ) {
        fortranio.SetFormat("6F12.7");
        result &= fortranio.ReadReal(Box.x);
        result &= fortranio.ReadReal(Box.y);
        result &= fortranio.ReadReal(Box.z);
        result &= fortranio.ReadReal(Box1.x);
        result &= fortranio.ReadReal(Box1.y);
        result &= fortranio.ReadReal(Box1.z);
    }
    if( (VelocitiesLoaded == false) && (Topology->BoxInfo.GetType() != AMBER_BOX_NONE) ) {
        if( loaded_v == 6 ) {
            Box.x = Velocities[0].x;
            Box.y = Velocities[0].y;
            Box.z = Velocities[0].z;
            Box1.x = Velocities[1].x;
            Box1.y = Velocities[1].y;
            Box1.z = Velocities[1].z;
            Velocities[0].x = 0.0;
            Velocities[0].y = 0.0;
            Velocities[0].z = 0.0;
            Velocities[1].x = 0.0;
            Velocities[1].y = 0.0;
            Velocities[1].z = 0.0;
        }
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberRestart::Save(FILE *fout)
{
    if( fout == NULL ) return(false);

    if( Topology == NULL ) {
        ES_ERROR("topology is not specified");
        return(false);
    }
    if( Topology->AtomList.GetNumberOfAtoms() == 0 ) {
        ES_ERROR("topology has no atom - is it loaded ?");
        return(false);
    }

    CFortranIO fortranio(fout);

    fortranio.SetFormat("20A4");

    // save title
    if( fortranio.WriteString(Title) == false ) {
        ES_ERROR("unable to save restart title");
        return(false);
    }
    fortranio.WriteEndOfSection();

    // fortranio.SetFormat("1I5");  // in AMBER6
    fortranio.SetFormat("1I6");    // in AMBER7
    if( fortranio.WriteInt(NumberOfAtoms) == false ) {
        ES_ERROR("unable to save number of atoms to restart file");
        return(false);
    }
    fortranio.ChangeFormat("1E15.7");
    if( fortranio.WriteReal(Time) == false ) {
        ES_ERROR("unable to save time position to restart file");
        return(false);
    }

    fortranio.WriteEndOfSection();

    fortranio.SetFormat("6F12.7");
    bool result = true;
    for(int i=0; i < NumberOfAtoms; i++) {
        result &= fortranio.WriteReal(Positions[i].x);
        result &= fortranio.WriteReal(Positions[i].y);
        result &= fortranio.WriteReal(Positions[i].z);
        if( result == false ) {
            CSmallString error;
            error << "unable to save position of atom " << i;
            ES_ERROR(error);
            return(false);
        }
    }
    fortranio.WriteEndOfSection();

    if( VelocitiesLoaded == true ) {
        fortranio.SetFormat("6F12.7");
        for(int i=0; i < NumberOfAtoms; i++) {
            result &= fortranio.WriteReal(Velocities[i].x);
            result &= fortranio.WriteReal(Velocities[i].y);
            result &= fortranio.WriteReal(Velocities[i].z);
            if( result == false ) {
                CSmallString error;
                error << "unable to save velocities of atom " << i;
                ES_ERROR(error);
                return(false);
            }
        }
        fortranio.WriteEndOfSection();
    }

    if( IsBoxPresent() == true ) {
        fortranio.SetFormat("6F12.7");
        result &= fortranio.WriteReal(Box.x);
        result &= fortranio.WriteReal(Box.y);
        result &= fortranio.WriteReal(Box.z);
        result &= fortranio.WriteReal(Box1.x);
        result &= fortranio.WriteReal(Box1.y);
        result &= fortranio.WriteReal(Box1.z);
        if( result == false ) {
            ES_ERROR("unable to save box info");
            return(false);
        }
        fortranio.WriteEndOfSection();
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberRestart::LoadSnapshot(CXMLElement* p_ele)
{
    if( p_ele == NULL ) {
        ES_ERROR("p_ele is NULL");
        return(false);
    }

    if( p_ele->GetName() != "SNAPSHOT" ) {
        ES_ERROR("p_ele is not SNAPSHOT");
        return(false);
    }

    if( Create() == false ) {
        ES_ERROR("unable to allocate memory for restart file");
        return(false);
    }

    if( Positions == NULL ) {
        ES_ERROR("Positions is NULL");
        return(false);
    }

    CXMLBinData* p_pele = p_ele->GetFirstChildBinData("POSITIONS");
    if( p_pele == NULL ) {
        ES_ERROR("unable to get BinData element POSITIONS");
        return(false);
    }

    // save box info
    if( Topology->BoxInfo.GetType() != AMBER_BOX_NONE ) {
        CXMLElement* p_bele = p_ele->GetFirstChildElement("BOX");
        if( p_bele == NULL ) {
            ES_ERROR("unable to get element BOX");
            return(false);
        }

        bool result = true;

        result &= p_bele->GetAttribute("x",Box.x);
        result &= p_bele->GetAttribute("y",Box.y);
        result &= p_bele->GetAttribute("z",Box.z);

        if( result == false ) {
            ES_ERROR("unable to load box attributes");
            return(false);
        }

        // box angles from topology
        Box1.x = Topology->BoxInfo.GetBoxBeta();
        Box1.y = Topology->BoxInfo.GetBoxBeta();
        Box1.z = Topology->BoxInfo.GetBoxBeta();
    }

    // get array
    double* p_array = (double*)p_pele->GetData();

    if( p_array == NULL ) {
        ES_ERROR("no data");
        return(false);
    }

    // check array size
    if( p_pele->GetLength() != 3*NumberOfAtoms*sizeof(double) ) {
        ES_ERROR("data length mismatch (are both systems the same?)");
        return(false);
    }

    // copy data
    int loc = 0;
    for(int i=0; i<NumberOfAtoms; i++) {
        Positions[i].x = p_array[loc++];
        Positions[i].y = p_array[loc++];
        Positions[i].z = p_array[loc++];
    }

    VelocitiesLoaded = true; // put zeros as velocities

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberRestart::LoadVelocities(CXMLElement* p_ele)
{
    if( p_ele == NULL ) {
        ES_ERROR("p_ele is NULL");
        return(false);
    }

    if( p_ele->GetName() != "SNAPSHOT" ) {
        ES_ERROR("p_ele is not SNAPSHOT");
        return(false);
    }

    if( Velocities == NULL ) {
        ES_ERROR("Velocities is NULL");
        return(false);
    }

    CXMLBinData* p_pele = p_ele->GetFirstChildBinData("POSITIONS");
    if( p_pele == NULL ) {
        ES_ERROR("unable to get BinData element POSITIONS");
        return(false);
    }

    // get array
    double* p_array = (double*)p_pele->GetData();

    if( p_array == NULL ) {
        ES_ERROR("no data");
        return(false);
    }

    // check array size
    if( p_pele->GetLength() != 3*NumberOfAtoms*sizeof(double) ) {
        ES_ERROR("data length mismatch (are both systems the same?)");
        return(false);
    }

    // copy data
    int loc = 0;
    for(int i=0; i<NumberOfAtoms; i++) {
        Velocities[i].x = p_array[loc++];
        Velocities[i].y = p_array[loc++];
        Velocities[i].z = p_array[loc++];
    }

    VelocitiesLoaded = true;

    return(true);
}

//------------------------------------------------------------------------------

void CAmberRestart::SaveSnapshot(CXMLElement* p_ele)
{
    if( p_ele == NULL ) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    if( p_ele->GetName() != "SNAPSHOT" ) {
        INVALID_ARGUMENT("p_ele is not SNAPSHOT");
    }

    if( Positions == NULL ) {
        LOGIC_ERROR("Positions are NULL")
    }

    CXMLBinData* p_pele = p_ele->CreateChildBinData("POSITIONS");

    // save box info
    if( Topology->BoxInfo.GetType() != AMBER_BOX_NONE ) {
        CXMLElement* p_bele = p_ele->CreateChildElement("BOX");

        p_bele->SetAttribute("x",Box.x);
        p_bele->SetAttribute("y",Box.y);
        p_bele->SetAttribute("z",Box.z);
    }

    p_pele->CopyData(&Positions[0].x,3*NumberOfAtoms*sizeof(double),EXBDT_DOUBLE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

const CPoint& CAmberRestart::GetPosition(int index) const
{
    if( Positions == NULL ) return(zero);
    return( Positions[index] );
}

//---------------------------------------------------------------------------

double CAmberRestart::GetMass(int index) const
{
    if( Topology == NULL ) return(0.0);
    return(Topology->AtomList.GetAtom(index)->GetMass());
}

//---------------------------------------------------------------------------

void CAmberRestart::SetPosition(int index,const CPoint& pos)
{
    if( Positions == NULL ) return;
    Positions[index] = pos;
}

//---------------------------------------------------------------------------

const CPoint& CAmberRestart::GetVelocity(int index) const
{
    if( Velocities == NULL ) return(zero);
    return( Velocities[index] );
}

//---------------------------------------------------------------------------

void CAmberRestart::SetVelocity(int index,const CPoint& vel)
{
    if( Velocities == NULL ) return;
    Velocities[index] = vel;
    VelocitiesLoaded = true;
}

//---------------------------------------------------------------------------

const CPoint& CAmberRestart::GetBox(void) const
{
    return(Box);
}

//---------------------------------------------------------------------------

const CPoint& CAmberRestart::GetAngles(void) const
{
    return(Box1);
}

//---------------------------------------------------------------------------

void CAmberRestart::SetBox(const CPoint& newbox)
{
    if( Topology->BoxInfo.GetType() == AMBER_BOX_NONE ) return;
    Box = newbox;
}

//---------------------------------------------------------------------------

void CAmberRestart::SetAngles(const CPoint& newangles)
{
    if( Topology->BoxInfo.GetType() == AMBER_BOX_NONE ) return;
    Box1 = newangles;
}

//---------------------------------------------------------------------------

bool CAmberRestart::AreVelocitiesLoaded(void) const
{
    return(VelocitiesLoaded);
}

//---------------------------------------------------------------------------

bool CAmberRestart::IsBoxPresent(void) const
{
    if( Topology->BoxInfo.GetType() == AMBER_BOX_NONE ) return(false);
    return(true);
}

//------------------------------------------------------------------------------

double CAmberRestart::GetTime(void)
{
    return(Time);
}

//------------------------------------------------------------------------------

void CAmberRestart::SetTime(const double time)
{
    Time = time;
}

//------------------------------------------------------------------------------

const CSmallString& CAmberRestart::GetTitle(void)
{
    return(Title);
}

//------------------------------------------------------------------------------

void CAmberRestart::SetTitle(const CSmallString& title)
{
    Title = title;
}

//------------------------------------------------------------------------------

void CAmberRestart::operator = (const CAmberRestart& src)
{
    if( GetNumberOfAtoms() != src.GetNumberOfAtoms() ){
        RUNTIME_ERROR("inconsistent number of atoms");
    }

    for(int i=0; i < GetNumberOfAtoms(); i++ ){
        CPoint pos = src.GetPosition(i);
        SetPosition(i,pos);
        if( src.AreVelocitiesLoaded() ){
            CPoint vel = src.GetVelocity(i);
            SetVelocity(i,vel);
        }
    }
    if( src.IsBoxPresent() ){
        SetBox(src.GetBox());
        SetAngles(src.GetAngles());
    }
}

//------------------------------------------------------------------------------

double* CAmberRestart::GetCoordinatesBuffer(void)
{
    if( Positions == NULL ) return(NULL);
    return( &(Positions->x) );
}

//------------------------------------------------------------------------------

double* CAmberRestart::GetVelocitiesBuffer(void)
{
    if( Velocities == NULL ) return(NULL);
    return( &(Velocities->x) );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
