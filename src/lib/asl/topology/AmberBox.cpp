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
#include <AmberBox.hpp>
#include <FortranIO.hpp>
#include <ErrorSystem.hpp>
#include <values.h>
#include <float.h>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberBox::CAmberBox(void)
{
    IFBOX = 0;
    NSP = NULL;
    IPTRES = 0;
    NSPM = 0;
    NSPSOL = 0;

    memset(DIMM,0,sizeof(DIMM));
    memset(ANGS,0,sizeof(ANGS));
    memset(UCELL,0,sizeof(UCELL));
    memset(RECIP,0,sizeof(RECIP));

    Volume = 0.0;
    Radius = 0.0;

    NumOfSoluteAtoms = 0;
}

//------------------------------------------------------------------------------

CAmberBox::~CAmberBox(void)
{
    FreeFields();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberBox::InitFields(int iIFBOX)
{
    IFBOX = iIFBOX;
    if( NSP != NULL ) delete[] NSP;
    NSP = NULL;
    IPTRES = 0;
    NSPM = 0;
    NSPSOL = 0;

    memset(DIMM,0,sizeof(DIMM));
    memset(ANGS,0,sizeof(ANGS));
    memset(UCELL,0,sizeof(UCELL));
    memset(RECIP,0,sizeof(RECIP));

    Volume = 0.0;
    Radius = 0.0;

    NumOfSoluteAtoms = 0;
}

//------------------------------------------------------------------------------

void CAmberBox::FreeFields(void)
{
    IFBOX = 0;
    if( NSP != NULL ) delete[] NSP;
    NSP = NULL;
    IPTRES = 0;
    NSPM = 0;
    NSPSOL = 0;

    memset(DIMM,0,sizeof(DIMM));
    memset(ANGS,0,sizeof(ANGS));
    memset(UCELL,0,sizeof(UCELL));
    memset(RECIP,0,sizeof(RECIP));

    NumOfSoluteAtoms = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

EAmberBoxType CAmberBox::GetType(void)
{
    switch(IFBOX) {
    case 0:
        return(AMBER_BOX_NONE);
    case 1:
        return(AMBER_BOX_STANDARD);
    case 2:
        return(AMBER_BOX_OCTAHEDRAL);
    default:
        return(AMBER_BOX_NONE);
    }
}

//------------------------------------------------------------------------------

void CAmberBox::SetType(EAmberBoxType type)
{
    switch(type) {
    case AMBER_BOX_NONE:
        IFBOX=0;
        break;
    case AMBER_BOX_STANDARD:
        IFBOX=1;
        break;
    case AMBER_BOX_OCTAHEDRAL:
        IFBOX=2;
        break;
    default:
        IFBOX=0;
        break;
    }
}

//------------------------------------------------------------------------------

int CAmberBox::GetNumberOfSoluteAtoms(void)
{
    int num=0;
    if( NSP == NULL ) return(num);

    for(int i=0; i<NSPSOL-1; i++) {
        num += NSP[i];
    }

    return(num);
}

//------------------------------------------------------------------------------

int CAmberBox::GetNumberOfSolventAtoms(void)
{
    int num=0;
    if( NSP == NULL ) return(num);

    for(int i=NSPSOL; i<NSPM-1; i++) {
        num += NSP[i];
    }

    return(num);
}

//------------------------------------------------------------------------------

int CAmberBox::GetNumberOfMolecules(void)
{
    return(NSPM);
}

//------------------------------------------------------------------------------

bool CAmberBox::SetNumberOfMolecules(int number)
{
    if( NSP != NULL ) delete[] NSP;
    NSP = NULL;
    IPTRES = 0;
    NSPM = 0;
    NSPSOL = 0;

    if( number > 0 ) {
        NSPM = number;
        NSP = new int[NSPM];
        if( NSP == NULL ) return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBox::SetFirstSolventMolecule(int index)
{
    if( (index < 0) || (index >=NSPM) ) return(false);
    NSPSOL = index + 1;
    return(true);
}

//------------------------------------------------------------------------------

int CAmberBox::GetNumberOfSoluteMolecules(void)
{
    return(NSPSOL-1);
}

//------------------------------------------------------------------------------

int CAmberBox::GetNumberOfSolventMolecules(void)
{
    return(NSPM-NSPSOL+1);
}

//------------------------------------------------------------------------------

int CAmberBox::GetNumberOfAtomsInMolecule(int mol_index)
{
    return(NSP[mol_index]);
}

//------------------------------------------------------------------------------

void CAmberBox::SetNumberOfAtomsInMolecule(int mol_index,int number)
{
    NSP[mol_index] = number;
}

//------------------------------------------------------------------------------

CPoint CAmberBox::GetBoxDimmensions(void)
{
    CPoint dim;

    dim.x = DIMM[0];
    dim.y = DIMM[1];
    dim.z = DIMM[2];

    return(dim);
}

//------------------------------------------------------------------------------

void CAmberBox::SetBoxDimmensions(const CPoint& dim)
{
    DIMM[0] = dim.x;
    DIMM[1] = dim.y;
    DIMM[2] = dim.z;
}

//------------------------------------------------------------------------------

double CAmberBox::GetBoxBeta(void)
{
    return(ANGS[0]);
}

//------------------------------------------------------------------------------

void CAmberBox::SetBoxBeta(double angle)
{
    ANGS[0] = angle;
    ANGS[1] = ANGS[0];
    ANGS[2] = ANGS[0];
}

//------------------------------------------------------------------------------

CPoint CAmberBox::GetBoxAngles(void)
{
    CPoint angs;

    angs.x = ANGS[0];
    angs.y = ANGS[1];
    angs.z = ANGS[2];

    return(angs);
}

//------------------------------------------------------------------------------

void CAmberBox::SetBoxAngles(const CPoint& angs)
{
    ANGS[0] = angs.x;
    ANGS[1] = angs.y;
    ANGS[2] = angs.z;
}

//------------------------------------------------------------------------------

void CAmberBox::UpdateBoxMatrices(void)
{
    memset(UCELL,0,sizeof(UCELL));
    memset(RECIP,0,sizeof(RECIP));

    Volume = 0.0;
    Radius = 0.0;

    double d2r = M_PI / 180.0;

    switch(GetType()) {
    default:
    case AMBER_BOX_NONE:
        // nothing to do
        return;

    case AMBER_BOX_STANDARD:
        UCELL[0][0] = DIMM[0];
        UCELL[1][1] = DIMM[1];
        UCELL[2][2] = DIMM[2];
        RECIP[0][0] = 1.0 / UCELL[0][0];
        RECIP[1][1] = 1.0 / UCELL[1][1];
        RECIP[2][2] = 1.0 / UCELL[2][2];

        Volume = UCELL[0][0]*UCELL[1][1]*UCELL[2][2];
        Radius = 0.5 * UCELL[0][0];
        if( (UCELL[1][1] < UCELL[2][2]) && (UCELL[1][1] < UCELL[0][0]) ) Radius = 0.5 * UCELL[1][1];
        if( (UCELL[2][2] < UCELL[1][1]) && (UCELL[2][2] < UCELL[0][0]) ) Radius = 0.5 * UCELL[2][2];
        break;

    case AMBER_BOX_OCTAHEDRAL: {
        double u12[3];
        double u20[3];
        double u01[3];

        // UCELL - orthogonalization matrix
        UCELL[0][0] = DIMM[0];
        UCELL[1][0] = 0.0;
        UCELL[2][0] = 0.0;

        UCELL[0][1] = DIMM[1]*cos(ANGS[2]*d2r);
        UCELL[1][1] = DIMM[1]*sin(ANGS[2]*d2r);
        UCELL[2][1] = 0.0;

        UCELL[0][2] = DIMM[2]*cos(ANGS[1]*d2r);
        UCELL[1][2] = (DIMM[1]*DIMM[2]*cos(ANGS[0]*d2r) - UCELL[0][2]*UCELL[0][1])/UCELL[1][1];
        UCELL[2][2] = sqrt(DIMM[2]*DIMM[2] - UCELL[0][2]*UCELL[0][2] - UCELL[1][2]*UCELL[1][2]);

        //         for(int i=0; i < 3; i++){
        //             for(int j=0; j < 3; j++){
        //                 printf("%f ",UCELL[i][j]);
        //                 }
        //             printf("\n");
        //             }

        // RECIP - deorthogonalization matrix
        u12[0] = UCELL[1][1]*UCELL[2][2] - UCELL[2][1]*UCELL[1][2];
        u12[1] = UCELL[2][1]*UCELL[0][2] - UCELL[0][1]*UCELL[2][2];
        u12[2] = UCELL[0][1]*UCELL[1][2] - UCELL[1][1]*UCELL[0][2];

        u20[0] = UCELL[1][2]*UCELL[2][0] - UCELL[2][2]*UCELL[1][0];
        u20[1] = UCELL[2][2]*UCELL[0][0] - UCELL[0][2]*UCELL[2][0];
        u20[2] = UCELL[0][2]*UCELL[1][0] - UCELL[1][2]*UCELL[0][0];

        u01[0] = UCELL[1][0]*UCELL[2][1] - UCELL[2][0]*UCELL[1][1];
        u01[1] = UCELL[2][0]*UCELL[0][1] - UCELL[0][0]*UCELL[2][1];
        u01[2] = UCELL[0][0]*UCELL[1][1] - UCELL[1][0]*UCELL[0][1];

        Volume = UCELL[0][0]*u12[0] + UCELL[1][0]*u12[1]+ UCELL[2][1]*u12[2];

        for(int i=0; i < 3; i++) RECIP[0][i] = u12[i] / Volume;
        for(int i=0; i < 3; i++) RECIP[1][i] = u20[i] / Volume;
        for(int i=0; i < 3; i++) RECIP[2][i] = u01[i] / Volume;

        //         for(int i=0; i < 3; i++){
        //             for(int j=0; j < 3; j++){
        //                 printf("%f ",RECIP[i][j]);
        //                 }
        //             printf("\n");
        //             }

        Radius = DIMM[0] * DIMM[1] * DIMM[2];

        for(int i=0; i < 3; i++) {
            double dist;
            dist = RECIP[i][0]*UCELL[0][i] + RECIP[i][1]*UCELL[1][i] + RECIP[i][2]*UCELL[2][i];
            dist = dist / sqrt( RECIP[i][0]*RECIP[i][0] + RECIP[i][1]*RECIP[i][1] + RECIP[i][2]*RECIP[i][2]);
            if( dist < Radius ) Radius = dist;
        }
        Radius = 0.5 * Radius;
    }
    break;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CAmberBox::GetVolume(void)
{
    return(Volume);
}

//------------------------------------------------------------------------------

double CAmberBox::GetLargestSphereRadius(void)
{
    return(Radius);
}

//------------------------------------------------------------------------------

const CPoint CAmberBox::GetBoxCenter(void)
{
    CPoint box_center;
    box_center.x = 0.5*UCELL[0][0] + 0.5*UCELL[0][1] + 0.5*UCELL[0][2];
    box_center.y = 0.5*UCELL[1][0] + 0.5*UCELL[1][1] + 0.5*UCELL[1][2];
    box_center.z = 0.5*UCELL[2][0] + 0.5*UCELL[2][1] + 0.5*UCELL[2][2];
    return(box_center);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

const CPoint CAmberBox::ImagePoint(const CPoint& pos,
        int kx,int ky,int kz,
        bool origin,bool familiar)
{
    CPoint ipos;

    // cell offset

    switch(GetType()) {
    default:
    case AMBER_BOX_NONE:
        return(pos);

    case AMBER_BOX_STANDARD: {
        double fx,fy,fz;
        // base cell offset
        fx = kx;
        fy = ky;
        fz = kz;

        // origin offset
        if( origin ) {
            fx -= 0.5;
            fy -= 0.5;
            fz -= 0.5;
        }
        ipos.x = pos.x - UCELL[0][0]*floor(pos.x*RECIP[0][0]) + fx*UCELL[0][0];
        ipos.y = pos.y - UCELL[1][1]*floor(pos.y*RECIP[1][1]) + fy*UCELL[1][1];
        ipos.z = pos.z - UCELL[2][2]*floor(pos.z*RECIP[2][2]) + fz*UCELL[2][2];
    }
    return(ipos);

    case AMBER_BOX_OCTAHEDRAL: {
        // calculate fractional coordinates
        double fx,fy,fz;
        fx = pos.x*RECIP[0][0] + pos.y*RECIP[0][1] + pos.z*RECIP[0][2];
        fy = pos.x*RECIP[1][0] + pos.y*RECIP[1][1] + pos.z*RECIP[1][2];
        fz = pos.x*RECIP[2][0] + pos.y*RECIP[2][1] + pos.z*RECIP[2][2];

        // in which are we?
        double ffx,ffy,ffz;
        ffx = floor(fx);
        ffy = floor(fy);
        ffz = floor(fz);

        // do familiar imaging
        if( familiar ) {
            double rfx,rfy,rfz;
            rfx = fx - floor(fx);
            rfy = fy - floor(fy);
            rfz = fz - floor(fz);

            // coordinates of cell centre
            double ox,oy,oz;
            ox = 0.5*UCELL[0][0] + 0.5*UCELL[0][1] + 0.5*UCELL[0][2];
            oy = 0.5*UCELL[1][0] + 0.5*UCELL[1][1] + 0.5*UCELL[1][2];
            oz = 0.5*UCELL[2][0] + 0.5*UCELL[2][1] + 0.5*UCELL[2][2];

            int lx,ly,lz;
            int mx=0,my=0,mz=0;

            // this is just upper guess of possible distance
            double min_dis = DBL_MAX;

            for(lx=-1; lx <= 1; lx++) {
                for(ly=-1; ly <= 1; ly++) {
                    for(lz=-1; lz <= 1; lz++) {
                        double ofx,ofy,ofz;
                        ofx = rfx + lx;
                        ofy = rfy + ly;
                        ofz = rfz + lz;

                        double px,py,pz;
                        px = ofx*UCELL[0][0] + ofy*UCELL[0][1] + ofz*UCELL[0][2];
                        py = ofx*UCELL[1][0] + ofy*UCELL[1][1] + ofz*UCELL[1][2];
                        pz = ofx*UCELL[2][0] + ofy*UCELL[2][1] + ofz*UCELL[2][2];

                        double ds = (px-ox)*(px-ox) + (py-oy)*(py-oy) + (pz-oz)*(pz-oz);

                        if( ds < min_dis ) {
                            mx = lx;
                            my = ly;
                            mz = lz;
                            min_dis = ds;
                        }
                    }
                }
            }
            ffx -= mx;
            ffy -= my;
            ffz -= mz;
        }

        // base cell offset
        ffx -= kx;
        ffy -= ky;
        ffz -= kz;

        // origin offset
        if( origin ) {
            ffx += 0.5;
            ffy += 0.5;
            ffz += 0.5;
        }

        double dx,dy,dz;

        dx = ffx*UCELL[0][0] + ffy*UCELL[0][1] + ffz*UCELL[0][2];
        dy = ffx*UCELL[1][0] + ffy*UCELL[1][1] + ffz*UCELL[1][2];
        dz = ffx*UCELL[2][0] + ffy*UCELL[2][1] + ffz*UCELL[2][2];
        // image point
        ipos.x = pos.x - dx;
        ipos.y = pos.y - dy;
        ipos.z = pos.z - dz;
    }
    return(ipos);
    }
}

//------------------------------------------------------------------------------

const CPoint CAmberBox::ImageVector(const CPoint& vec)
{
    CPoint ipos = GetBoxCenter() + vec;
    ipos = ImagePoint(ipos,0,0,0,false,true);
    ipos = ipos - GetBoxCenter();
    return(ipos);
}

//------------------------------------------------------------------------------

void CAmberBox::operator = (const CAmberBox& src)
{
//    void InitFields(int iIFBOX);

    InitFields(src.IFBOX);

    IPTRES = src.IPTRES;
    NSPM = src.NSPM;
    NSPSOL = src.NSPSOL;

    if(NSPM > 0) {
        NSP = new int[NSPM];
    }

    for(int i=0; i<NSPM; i++) {
        NSP[i] = src.NSP[i];
    }

    for(int i=0; i < 3;i++){
        DIMM[i] = src.DIMM[i];
        ANGS[i] = src.ANGS[i];
    }

    for(int i=0; i < 3;i++){
        for(int j=0; j < 3;j++){
            UCELL[i][j] = src.UCELL[i][j];
            RECIP[i][j] = src.RECIP[i][j];
        }
    }

    Volume = src.Volume;
    Radius = src.Radius;
    NumOfSoluteAtoms = src.NumOfSoluteAtoms;

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberBox::JoinSoluteMolecules(int number_of_first_molecules)
{
    if( (NSPM == 0)||(NSP == NULL) ) {
        ES_ERROR("box info is not present");
        return(false);
    }

    for(int i=1; i<number_of_first_molecules; i++) {
        NSP[0] += NSP[i];
    }

    // move data down
    for(int i=number_of_first_molecules; i<NSPM; i++) {
        NSP[i-number_of_first_molecules+1] = NSP[i];
    }

    // update indexes
    NSPM   -= (number_of_first_molecules-1);
    NSPSOL -= (number_of_first_molecules-1);

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberBox::LoadSolventPointers(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    if( fortranio.ReadInt(IPTRES) == false ) {
        ES_ERROR("unable load IPTRES item");
        return(false);
    }

    if( fortranio.ReadInt(NSPM) == false ) {
        ES_ERROR("unable load NSPM item");
        return(false);
    }

    if( fortranio.ReadInt(NSPSOL) == false ) {
        ES_ERROR("unable load NSPSOL item");
        return(false);
    }

    // allocate fields
    if( (NSP = new int[NSPM]) == NULL ) {
        ES_ERROR("unable allocate NSP array");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBox::LoadNumsOfMolecules(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    int* p_item = NSP;

    for(int i=0; i<NSPM; i++) {
        if( fortranio.ReadInt(*p_item) == false ) {
            ES_ERROR("unable load NSP item");
            return(false);
        }
        p_item++;
    }

    NumOfSoluteAtoms=0;
    for(int i=0; i<NSPSOL-1; i++) {
        NumOfSoluteAtoms += NSP[i];
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBox::LoadBoxInfo(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    if( fortranio.ReadReal(ANGS[0]) == false ) {
        ES_ERROR("unable load BETA item");
        return(false);
    }

    ANGS[1] = ANGS[0];
    ANGS[2] = ANGS[0];

    for(int i=0; i<3; i++) {
        if( fortranio.ReadReal(DIMM[i]) == false ) {
            ES_ERROR("nable load BOX item");
            return(false);
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBox::SaveAmber6Box(FILE* p_file)
{
    if( SaveSolventPointers(p_file,"12I6") == false ) return(false);
    if( SaveNumsOfMolecules(p_file,"12I6") == false ) return(false);
    if( SaveBoxInfo(p_file,"5E16.8") == false ) return(false);
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBox::SaveSolventPointers(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    if( fortranio.WriteInt(IPTRES) == false ) {
        ES_ERROR("unable save IPTRES item");
        return(false);
    }

    if( fortranio.WriteInt(NSPM) == false ) {
        ES_ERROR("unable save NSPM item");
        return(false);
    }

    if( fortranio.WriteInt(NSPSOL) == false ) {
        ES_ERROR("unable save NSPSOL item");
        return(false);
    }

    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBox::SaveNumsOfMolecules(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    int* p_item = NSP;

    for(int i=0; i<NSPM; i++) {
        if( fortranio.WriteInt(*p_item) == false ) {
            ES_ERROR("unable save NSP item");
            return(false);
        }
        p_item++;
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//------------------------------------------------------------------------------

bool CAmberBox::SaveBoxInfo(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    if( fortranio.WriteReal(ANGS[0]) == false ) {
        ES_ERROR("unable save BETA item");
        return(false);
    }

    for(int i=0; i<3; i++) {
        if( fortranio.WriteReal(DIMM[i]) == false ) {
            ES_ERROR("unable save BOX item");
            return(false);
        }
    }
    fortranio.WriteEndOfSection();
    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

