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
#include <errno.h>
#include <AmberTopology.hpp>
#include <FortranIO.hpp>
#include <ErrorSystem.hpp>
#include <list>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberTopology::CAmberTopology(void)
    : ResidueList(this)
{
    FakeTopology = false;
    Version = AMBER_VERSION_7;
    NHPARM = 0;
    NPARM = 0;

    // amber_7 defaults
    fTITLE="20a4";
    fPOINTERS="10I8";
    fATOM_NAME="20a4";
    fCHARGE="5E16.8";
    fMASS="5E16.8";
    fATOM_TYPE_INDEX="10I8";
    fNUMBER_EXCLUDED_ATOMS="10I8";
    fNONBONDED_PARM_INDEX="10I8";
    fRESIDUE_LABEL="20a4";
    fRESIDUE_POINTER="10I8";
    fBOND_FORCE_CONSTANT="5E16.8";
    fBOND_EQUIL_VALUE="5E16.8";
    fANGLE_FORCE_CONSTANT="5E16.8";
    fANGLE_EQUIL_VALUE="5E16.8";
    fDIHEDRAL_FORCE_CONSTANT="5E16.8";
    fDIHEDRAL_PERIODICITY="5E16.8";
    fDIHEDRAL_PHASE="5E16.8";
    fSOLTY="5E16.8";
    fLENNARD_JONES_ACOEF="5E16.8";
    fLENNARD_JONES_BCOEF="5E16.8";
    fBONDS_INC_HYDROGEN="10I8";
    fBONDS_WITHOUT_HYDROGEN="10I8";
    fANGLES_INC_HYDROGEN="10I8";
    fANGLES_WITHOUT_HYDROGEN="10I8";
    fDIHEDRALS_INC_HYDROGEN="10I8";
    fDIHEDRALS_WITHOUT_HYDROGEN="10I8";
    fEXCLUDED_ATOMS_LIST="10I8";
    fHBOND_ACOEF="5E16.8";
    fHBOND_BCOEF="5E16.8";
    fHBCUT="5E16.8";
    fAMBER_ATOM_TYPE="20a4";
    fTREE_CHAIN_CLASSIFICATION="20a4";
    fJOIN_ARRAY="10I8";
    fIROTAT="10I8";
    fSOLVENT_POINTERS="10I8";
    fATOMS_PER_MOLECULE="10I8";
    fBOX_DIMENSIONS="5E16.8";
    fRADIUS_SET="1A80";
    fRADII="5E16.8";
    fSCREEN="5E16.8";
    fPERT_BOND_ATOMS="";
    fPERT_BOND_PARAMS="";
    fPERT_ANGLE_ATOMS="";
    fPERT_ANGLE_PARAMS="";
    fPERT_DIHEDRAL_ATOMS="";
    fPERT_DIHEDRAL_PARAMS="";
    fPERT_RESIDUE_NAME="";
    fPERT_ATOM_NAME="";
    fPERT_ATOM_SYMBOL="";
    fALMPER="";
    fIAPER="";
    fPERT_ATOM_TYPE_INDEX="";
    fPERT_CHARGE="";
    fSCEE_SCALE_FACTOR = "5E16.8";
    fSCNB_SCALE_FACTOR = "5E16.8";
    fATOMIC_NUMBER = "10I8";
}

//---------------------------------------------------------------------------

CAmberTopology::~CAmberTopology(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberTopology::Clean(void)
{
    Version = AMBER_VERSION_NONE;
    NHPARM = 0;
    NPARM = 0;
    Name = NULL;

    AtomList.FreeFields();
    ResidueList.FreeFields();
    BondList.FreeFields();
    AngleList.FreeFields();
    DihedralList.FreeFields();
    NonBondedList.FreeFields();
    BoxInfo.FreeFields();
    CapInfo.FreeFields();

    FakeTopology = false;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTopology::InitResidueNumOfBondsWithHydrogen(void)
{
    /* very slow implementation
    for(int j=0;j<BondList.GetNumberOfBondsWithHydrogen();j++){
       // find residue
       int at = BondList.GetBondWithHydrogen(j)->GetIB();
       for(int i = 0; i < ResidueList.GetNumberOfResidues(); i++){
           CAmberResidue* p_res = ResidueList.GetResidue(i);
           int pos = at - p_res->GetFirstAtomIndex();
           if( (pos >= 0) && (p_res->GetNumberOfAtoms() > pos) ){
               p_res->NumOfBondsWithHydrogen++;
               break;
               }
           p_res++;
           }
       }
    */
    /*
    int bond_id = 0;


    // we assume that bond list is sorted !!!!
    for(int i = 0; i < ResidueList.GetNumberOfResidues(); i++){
       CAmberResidue* p_res = ResidueList.GetResidue(i);
       int start = p_res->GetFirstAtomIndex();
       int stop = start + p_res->GetNumberOfAtoms();
       bool passed = bond_id < BondList.GetNumberOfBondsWithHydrogen();
       while( passed == true ){
           CAmberBond* p_bond = BondList.GetBondWithHydrogen(bond_id);
           if( (p_bond->GetIB() >= start) && (p_bond->GetIB() < stop) ){
               p_res->NumOfBondsWithHydrogen++;
               }
           else{
               if( (p_bond->GetJB() >= start) && (p_bond->GetJB() < stop) ){
                   p_res->NumOfBondsWithHydrogen++;
                   }
                 else{
                     passed = false;
                   break;
                   }
               }
           bond_id++;
           if( bond_id >= BondList.GetNumberOfBondsWithHydrogen() ) passed = false;
           }
       }
    */

    for(int j=0; j<BondList.GetNumberOfBondsWithHydrogen(); j++) {
        // find residue
        int at1 = BondList.GetBondWithHydrogen(j)->GetIB();
        int at2 = BondList.GetBondWithHydrogen(j)->GetJB();
        CAmberAtom* p_at1 = AtomList.GetAtom(at1);
        CAmberAtom* p_at2 = AtomList.GetAtom(at2);
        if( (p_at1 != NULL) && (p_at2 != NULL) ) {
            if( p_at1->GetResidue() == p_at2->GetResidue() ) {
                if( p_at1->GetResidue() != NULL ) p_at1->GetResidue()->NumOfBondsWithHydrogen++;
            } else {
                ES_ERROR("bond with hydrogen between two residues");
                return(false);
            }
        }
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberTopology::InitMoleculeIndexes(void)
{
    int            mol_index = 0;
    CAmberResidue* p_last_res = NULL;

    // assign initial indexes, we assume that residues are not built from several molecules

    for(int j=0; j<AtomList.GetNumberOfAtoms(); j++) {
        CAmberAtom* p_atom = AtomList.GetAtom(j);
        if( p_atom != NULL ) {
            if( p_last_res != p_atom->GetResidue() ) {
                p_last_res = p_atom->GetResidue();
                mol_index++;
            }
            p_atom->MoleculeIndex = mol_index;
        }
    }

    CAmberAtom*     p_at1;
    CAmberAtom*     p_at2;
    bool            passed = false;
    int             pass_index = 0;

    // this is iterative aproach, but it is quit fast

    while( passed == false ) {
        // printf("Pass index: %d\n",pass_index);
        pass_index++;
        passed = true;

        for(int i = 0; i < BondList.GetNumberOfBondsWithoutHydrogen(); i++) {
            CAmberBond* p_bond = BondList.GetBondWithoutHydrogen(i);
            if( p_bond != NULL ) {

                p_at1 = AtomList.GetAtom(p_bond->GetIB());
                p_at2 = AtomList.GetAtom(p_bond->GetJB());

                if( p_at1->GetMoleculeIndex() != p_at2->GetMoleculeIndex() ) {
                    passed = false;
                    // use lower index
                    if( p_at1->GetMoleculeIndex() < p_at2->GetMoleculeIndex() ) {
                        p_at2->MoleculeIndex = p_at1->GetMoleculeIndex();
                    } else {
                        p_at1->MoleculeIndex = p_at2->GetMoleculeIndex();
                    }
                }

            }
        }
    }

    // check consistency
    for(int i = 0; i < BondList.GetNumberOfBondsWithHydrogen(); i++) {
        CAmberBond* p_bond = BondList.GetBondWithHydrogen(i);
        if( p_bond != NULL ) {

            p_at1 = AtomList.GetAtom(p_bond->GetIB());
            p_at2 = AtomList.GetAtom(p_bond->GetJB());

            if( p_at1->GetMoleculeIndex() != p_at2->GetMoleculeIndex() ) {
                // use lower index
                if( p_at1->GetMoleculeIndex() < p_at2->GetMoleculeIndex() ) {
                    p_at2->MoleculeIndex = p_at1->GetMoleculeIndex();
                } else {
                    p_at1->MoleculeIndex = p_at2->GetMoleculeIndex();
                }
            }

        }
    }

    return(true);
}

//------------------------------------------------------------------------------

void CAmberTopology::BuidListOfNeighbourAtoms(void)
{
    for(int i = 0; i < BondList.GetNumberOfBonds(); i++) {
        CAmberBond* p_bond = BondList.GetBond(i);
        int at1 = p_bond->GetIB();
        int at2 = p_bond->GetJB();

        CAmberAtom* p_at1 = AtomList.GetAtom(at1);
        p_at1->NeighbourIndexes.insert(at2);

        CAmberAtom* p_at2 = AtomList.GetAtom(at2);
        p_at2->NeighbourIndexes.insert(at1);
    }
}

//------------------------------------------------------------------------------

void CAmberTopology::operator = (const CAmberTopology& src)
{
    Clean();

    /// currently loaded version of topology
    Version = src.Version;

    /// name of topology
    Name = src.Name;

    /// ITITL  : title
    ITITL = src.ITITL;

    /// NHPARM : currently not used
    NHPARM = src.NHPARM;

    /// NPARM  : currently not used
    NPARM = src.NPARM;

    // local copy of formats
    fTITLE = src.fTITLE;
    fPOINTERS = src.fPOINTERS;
    fATOM_NAME = src.fATOM_NAME;
    fCHARGE = src.fCHARGE;
    fMASS = src.fMASS;
    fATOM_TYPE_INDEX = src.fATOM_TYPE_INDEX;
    fNUMBER_EXCLUDED_ATOMS = src.fNUMBER_EXCLUDED_ATOMS;
    fNONBONDED_PARM_INDEX = src.fNONBONDED_PARM_INDEX;
    fRESIDUE_LABEL = src.fRESIDUE_LABEL;
    fRESIDUE_POINTER = src.fRESIDUE_POINTER;
    fBOND_FORCE_CONSTANT = src.fBOND_FORCE_CONSTANT;
    fBOND_EQUIL_VALUE = src.fBOND_EQUIL_VALUE;
    fANGLE_FORCE_CONSTANT = src.fANGLE_FORCE_CONSTANT;
    fANGLE_EQUIL_VALUE = src.fANGLE_EQUIL_VALUE;
    fDIHEDRAL_FORCE_CONSTANT = src.fDIHEDRAL_FORCE_CONSTANT;
    fDIHEDRAL_PERIODICITY = src.fDIHEDRAL_PERIODICITY;
    fDIHEDRAL_PHASE = src.fDIHEDRAL_PHASE;
    fSOLTY = src.fSOLTY;
    fLENNARD_JONES_ACOEF = src.fLENNARD_JONES_ACOEF;
    fLENNARD_JONES_BCOEF = src.fLENNARD_JONES_BCOEF;
    fBONDS_INC_HYDROGEN = src.fBONDS_INC_HYDROGEN;
    fBONDS_WITHOUT_HYDROGEN = src.fBONDS_WITHOUT_HYDROGEN;
    fANGLES_INC_HYDROGEN = src.fANGLES_INC_HYDROGEN;
    fANGLES_WITHOUT_HYDROGEN = src.fANGLES_WITHOUT_HYDROGEN;
    fDIHEDRALS_INC_HYDROGEN = src.fDIHEDRALS_INC_HYDROGEN;
    fDIHEDRALS_WITHOUT_HYDROGEN = src.fDIHEDRALS_WITHOUT_HYDROGEN;
    fEXCLUDED_ATOMS_LIST = src.fEXCLUDED_ATOMS_LIST;
    fHBOND_ACOEF = src.fHBOND_ACOEF;
    fHBOND_BCOEF = src.fHBOND_BCOEF;
    fHBCUT = src.fHBCUT;
    fAMBER_ATOM_TYPE = src.fAMBER_ATOM_TYPE;
    fTREE_CHAIN_CLASSIFICATION = src.fTREE_CHAIN_CLASSIFICATION;
    fJOIN_ARRAY = src.fJOIN_ARRAY;
    fIROTAT = src.fIROTAT;
    fSOLVENT_POINTERS = src.fSOLVENT_POINTERS;
    fATOMS_PER_MOLECULE = src.fATOMS_PER_MOLECULE;
    fBOX_DIMENSIONS = src.fBOX_DIMENSIONS;
    fRADIUS_SET = src.fRADIUS_SET;
    fRADII = src.fRADII;
    fSCREEN = src.fSCREEN;
    fPERT_BOND_ATOMS = src.fPERT_BOND_ATOMS;
    fPERT_BOND_PARAMS = src.fPERT_BOND_PARAMS;
    fPERT_ANGLE_ATOMS = src.fPERT_ANGLE_ATOMS;
    fPERT_ANGLE_PARAMS = src.fPERT_ANGLE_PARAMS;
    fPERT_DIHEDRAL_ATOMS = src.fPERT_DIHEDRAL_ATOMS;
    fPERT_DIHEDRAL_PARAMS = src.fPERT_DIHEDRAL_PARAMS;
    fPERT_RESIDUE_NAME = src.fPERT_RESIDUE_NAME;
    fPERT_ATOM_NAME = src.fPERT_ATOM_NAME;
    fPERT_ATOM_SYMBOL = src.fPERT_ATOM_SYMBOL;
    fALMPER = src.fALMPER;
    fIAPER = src.fIAPER;
    fPERT_ATOM_TYPE_INDEX = src.fPERT_ATOM_TYPE_INDEX;
    fPERT_CHARGE = src.fPERT_CHARGE;

    fSCEE_SCALE_FACTOR = src.fSCEE_SCALE_FACTOR;
    fSCNB_SCALE_FACTOR = src.fSCNB_SCALE_FACTOR;
    fATOMIC_NUMBER = src.fATOMIC_NUMBER;

    AtomList = src.AtomList;
    ResidueList = src.ResidueList;
    BondList = src.BondList;
    AngleList = src.AngleList;
    DihedralList = src.DihedralList;
    NonBondedList = src.NonBondedList;
    BoxInfo = src.BoxInfo;
    CapInfo = src.CapInfo;

    FakeTopology = src.FakeTopology;

    ResidueList.ReinitAtomResiduePointers(&AtomList);
    InitMoleculeIndexes();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberTopology::PrintInfo(bool short_info,FILE* p_out)
{
    if( p_out == NULL ) {
        p_out = stdout;
    }

    BoxInfo.UpdateBoxMatrices();

    fprintf(p_out,"\n");

    if( Name != NULL ){
        fprintf(p_out," Topology name       : %s\n",(const char*)Name);
    } else {
        fprintf(p_out," Topology name       : -none-\n");
    }
    fprintf(p_out," Topology title      : %s\n",(const char*)ITITL);
    fprintf(p_out," Format version      : ");
    switch(Version) {
    case AMBER_VERSION_NONE:
        fprintf(p_out,"no version is provided\n");
        break;
    case AMBER_VERSION_6:
        fprintf(p_out,"Amber 6\n");
        break;
    case AMBER_VERSION_7:
        fprintf(p_out,"Amber 7\n");
        break;
    }

    if( short_info == true ) {
        if( FakeTopology ){
        fprintf(p_out," Completness         : -fake-\n");
        } else {
        fprintf(p_out," Completness         : -full-\n");
        }
        fprintf(p_out," Number of atoms     : % 8d\n",AtomList.NATOM);
        fprintf(p_out," Number of residues  : % 8d\n",ResidueList.NRES);
        if( AtomList.RadiusSet != NULL ) {
            fprintf(p_out," Radius set          : %s\n",(const char*)AtomList.RadiusSet);
        } else {
            fprintf(p_out," Radius set          : unknown/not present\n");
        }

        if( BoxInfo.GetType() != AMBER_BOX_NONE ) {
            switch( BoxInfo.GetType() ) {
            case AMBER_BOX_STANDARD:
                fprintf(p_out," System box            : rectangular\n");
                break;
            case AMBER_BOX_OCTAHEDRAL:
                fprintf(p_out," System box          : truncated octahedral\n");
                break;
            case AMBER_BOX_NONE:
                fprintf(p_out," System box          : none\n");
                break;
            }
            fprintf(p_out,"\n");
        } else {
            fprintf(p_out," System box          : none\n");
            fprintf(p_out,"\n");
        }
        return;
    }

    fprintf(p_out,"================================================================================\n");

    // long form info -------------------------------------------------------
    if( FakeTopology ){
    fprintf(p_out," Completness         : -fake-\n");
    } else {
    fprintf(p_out," Completness         : -full-\n");
    }
    fprintf(p_out," Number of atoms     : % 8d\n",AtomList.NATOM);
    fprintf(p_out," Number of bonds     : % 8d\n",BondList.NBONH+BondList.MBONA);
    fprintf(p_out," Number of X-H bonds : % 8d\n",BondList.NBONH);
    fprintf(p_out," Number of angles    : % 8d\n",AngleList.NTHETH+AngleList.MTHETA);
    fprintf(p_out," Number of dihedrals : % 8d\n",DihedralList.NPHIH+DihedralList.MPHIA);
    if( DihedralList.SCEEFactorsLoaded ){
    fprintf(p_out," SCEE factors loaded : yes\n");
    } else {
    fprintf(p_out," SCEE factors loaded : no\n");
    }
    if( DihedralList.SCNBFactorsLoaded ){
    fprintf(p_out," SCNB factors loaded : yes\n");
    } else {
    fprintf(p_out," SCNB factors loaded : no\n");
    }
    fprintf(p_out," Number of freedom   : % 8d  (with SHAKE : % 7d)\n",AtomList.NATOM*3-6,AtomList.NATOM*3-6-BondList.NBONH);
    fprintf(p_out," Number of residues  : % 8d\n",ResidueList.NRES);

    double total_charge = 0.0;
    for(int i = 0; i < AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = AtomList.GetAtom(i);
        if( p_atom != NULL ) {
            total_charge += p_atom->GetStandardCharge();
        }
    }
    if( AtomList.RadiusSet != NULL ) {
        fprintf(p_out," Radius set          : %s\n",(const char*)AtomList.RadiusSet);
    } else {
        fprintf(p_out," Radius set          : unknown/not present\n");
    }
    fprintf(p_out," Total charge        : % 8.4f\n",total_charge);


    if( AtomList.IFPERT == 1 ) {
        int pertatoms = 0;
        for(int i = 0; i < AtomList.GetNumberOfAtoms(); i++) {
            CAmberAtom* p_atom = AtomList.GetAtom(i);
            if( p_atom != NULL ) {
                if( p_atom->IsPerturbed() == true ) pertatoms++;
            }
        }
        fprintf(p_out,"\n");
        fprintf(p_out," System contains perturbation info.\n");
        fprintf(p_out," -------------------------------------------\n");
        fprintf(p_out," Number of perturbed atoms     : % 8d\n",AtomList.NATOM);
        fprintf(p_out," Number of perturbed bonds     : % 8d\n",BondList.NBPER);
        fprintf(p_out," Number of perturbed angles    : % 8d\n",AngleList.MGPER);
        fprintf(p_out," Number of perturbed dihedrals : % 8d\n",DihedralList.NDPER);
    }
    fprintf(p_out,"\n");
    fprintf(p_out,"================================================================================\n");

    if( AtomList.HasPertInfo() == true ) {
        fprintf(p_out,"\n");
        fprintf(p_out," System contains perturbation information. \n");
        fprintf(p_out,"\n");
        fprintf(p_out,"================================================================================\n");
    }

    if( BoxInfo.GetType() != AMBER_BOX_NONE ) {
        fprintf(p_out,"\n");
        switch( BoxInfo.GetType() ) {
        case AMBER_BOX_STANDARD:
            fprintf(p_out," Rectangular box is present.\n");
            break;
        case AMBER_BOX_OCTAHEDRAL:
            fprintf(p_out," Truncated octahedral box is present.\n");
            break;
        case AMBER_BOX_NONE:
            break;
        }
        fprintf(p_out," -------------------------------------------\n");
        fprintf(p_out," a                         = %12.4f A\n",BoxInfo.DIMM[0]);
        fprintf(p_out," b                         = %12.4f A\n",BoxInfo.DIMM[1]);
        fprintf(p_out," c                         = %12.4f A\n",BoxInfo.DIMM[2]);
        fprintf(p_out," beta                      = %12.4f deg\n",BoxInfo.ANGS[0]);
        fprintf(p_out," Box volume                : %12.4f A^3\n",BoxInfo.GetVolume());
        fprintf(p_out," Largest sphere radius     : %12.4f A\n",BoxInfo.GetLargestSphereRadius());
        fprintf(p_out," Number of molecules       : % 12d\n",BoxInfo.NSPM);
        fprintf(p_out," First solvent molecule    : % 12d\n",BoxInfo.NSPSOL);
        fprintf(p_out," Number of atoms in each solute molecule:\n");

        CFortranIO fortranio(p_out);
        fortranio.SetFormat("12I5");

        for(int i=0; i<BoxInfo.NSPSOL-1; i++) {
            fortranio.WriteInt(BoxInfo.NSP[i]);
        }
        fortranio.WriteEndOfSection();
    } else {
        fprintf(p_out,"\n");
        fprintf(p_out," Box is NOT present. \n");
    }

    fprintf(p_out,"\n");
    fprintf(p_out,"================================================================================\n");

    fprintf(p_out,"\n");
    fprintf(p_out," Residue summary\n");
    fprintf(p_out,"--------------------------------------------------------------------------------\n");

    ResidueList.PrepareSortedResidues();

    CAmberResidue*  p_prev = NULL;
    CAmberResidue*  p_res = NULL;
    CAmberResidue*  p_next = NULL;
    char            buffer[15];
    int             rescnt = 0;

    CFortranIO fortranio(p_out);
    fortranio.SetFormat("4a15");

    for(int i = 0; i < ResidueList.GetNumberOfResidues(); i++) {
        p_res = ResidueList.GetSortedResidue(i);
        if( (p_prev == NULL) || (strncmp(p_prev->GetName(),p_res->GetName(),4) != 0 ) ) {
            if( p_prev != NULL ) {
                sprintf(buffer,"   %-4s[%05d]",p_prev->GetName(),rescnt);
                fortranio.WriteString(buffer);
                rescnt = 0;
            }
            p_prev = p_res;
        }
        rescnt++;
    }
    if( p_prev != NULL ) {
        sprintf(buffer,"   %-4s[%05d]",p_prev->GetName(),rescnt);
        fortranio.WriteString(buffer);
    }
    fortranio.WriteEndOfSection();

    ResidueList.FreeSortedResidueList();

    fprintf(p_out,"\n");
    fprintf(p_out,"================================================================================\n");
    fprintf(p_out,"\n");
    fprintf(p_out," List of Residues\n");
    fprintf(p_out,"--------------------------------------------------------------------------------\n");

    p_prev = NULL;
    p_res = NULL;
    p_next = NULL;
    rescnt = 0;

    fortranio.SetFormat("5a12");

    if( ResidueList.GetNumberOfResidues() != 0 ) p_res = ResidueList.GetResidue(0);
    if( p_res != NULL ) {
        sprintf(buffer," %5d %4s",p_res->GetIndex()+1,p_res->GetName());
        fortranio.WriteString(buffer);
        p_prev = p_res;
    }
    p_res = NULL;
    for(int i = 1; i < ResidueList.GetNumberOfResidues()-1; i++) {
        p_res = ResidueList.GetResidue(i);
        if( strcmp(p_res->GetName(),p_prev->GetName()) != 0 ) {
            sprintf(buffer," %5d %4s",p_res->GetIndex()+1,p_res->GetName());
            fortranio.WriteString(buffer);
            rescnt = 0;
            p_prev = p_res;
        } else {
            rescnt++;
            p_next = ResidueList.GetResidue(i+1);
            if( strcmp(p_res->GetName(),p_next->GetName()) != 0 ) {
                if( rescnt >= 3 ) {
                    sprintf(buffer," <= %4s=> ",p_res->GetName());
                    fortranio.WriteString(buffer);
                }
                sprintf(buffer," %5d %4s",p_res->GetIndex()+1,p_res->GetName());
                fortranio.WriteString(buffer);
                rescnt = 0;
                p_prev = p_res;
            } else {
                rescnt++;
            }
        }
    }
    if( ResidueList.GetNumberOfResidues() >= 2 ) p_res = ResidueList.GetResidue(ResidueList.GetNumberOfResidues()-1);
    if( p_res != NULL ) {
        if( strcmp(p_res->GetName(),p_prev->GetName()) == 0 ) {
            if( rescnt >= 3 ) {
                sprintf(buffer," <= %4s=> ",p_res->GetName());
                fortranio.WriteString(buffer);
            }
        }
        sprintf(buffer," %5d %4s",p_res->GetIndex()+1,p_res->GetName());
        fortranio.WriteString(buffer);
        p_prev = p_res;
    }
    fortranio.WriteEndOfSection();
    fprintf(p_out,"\n");
}

//------------------------------------------------------------------------------

void CAmberTopology::PrintAtoms(FILE* p_out)
{
    if( p_out == NULL ) {
        p_out = stdout;
    }

    fprintf(p_out,"\n");
    fprintf(p_out,"# SerIndex Name SeqIdx ResName Type  Charge \n");
    fprintf(p_out,"# -------- ---- ------ ------- ---- --------\n");

    for(int i=0; i < AtomList.GetNumberOfAtoms(); i++ ){
        CAmberAtom* p_atom = AtomList.GetAtom(i);
        fprintf(p_out,"%10d %4s %6d %4s    %4s %8.4f\n",p_atom->GetAtomIndex()+1,p_atom->GetName(),
                                                      p_atom->GetResidue()->GetIndex()+1,p_atom->GetResidue()->GetName(),
                                                      p_atom->GetType(),p_atom->GetStandardCharge());
    }
    fprintf(p_out,"\n");

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

EAmberVersion CAmberTopology::GetVersion(void)
{
    return(Version);
}

//---------------------------------------------------------------------------

CSmallString CAmberTopology::GetTitle(void)
{
    return(ITITL);
}

//---------------------------------------------------------------------------

void CAmberTopology::SetTitle(const CSmallString& title)
{
    ITITL = title;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTopology::Load(const CSmallString& file_name,bool allow_stdin)
{
    FILE* p_top;

    if( (allow_stdin == true) && (file_name == "-") ) {
        p_top = stdin;
    } else {
        p_top = fopen(file_name,"r");
        if( p_top == NULL ) {
            CSmallString error;
            error << "unable to open topology file '" << file_name << "' ("
                  << strerror(errno) << ")";
            ES_ERROR(error);
            return(false);
        }
    }

    bool result = Load(p_top);

    if( result == true ) {
        Name = file_name;
    }

    if( ! ((allow_stdin == true) && (file_name == "-")) ) {
        fclose(p_top);
    }

    return(result);
}

//------------------------------------------------------------------------------

bool CAmberTopology::Load(FILE* p_fin)
{
    Clean();

    char buffer[255];
    fgets(buffer,254,p_fin);

    bool result;

    if( strstr(buffer,"%VERSION") != NULL ) {
        Version=AMBER_VERSION_7;
        result = LoadAmber7(p_fin);
    } else {
        buffer[254]='\0';
        int pos = strlen(buffer) - 1;
        if( pos < 0 ) pos = 0;
        buffer[pos]='\0'; // remove \n character
        ITITL = buffer;
        Version=AMBER_VERSION_6;
        result = LoadAmber6(p_fin);
    }

    if( result == false ) {
        switch(Version) {
        case AMBER_VERSION_6:
            ES_ERROR("unable to load topology file in AMBER 6 format");
            break;
        case AMBER_VERSION_7:
        default:
            ES_ERROR("unable to load topology file in AMBER 7 format");
            break;
        }
    }

    return(result);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

class CFakeResidue{
public:
    CSmallString        Name;
    list<CSmallString>  Atoms;
};

//------------------------------------------------------------------------------

bool CAmberTopology::LoadFakeTopologyFromPDB(const CSmallString& file_name,bool mangle_names,bool allow_stdin)
{
    FILE* p_top;

    if( (allow_stdin == true) && (file_name == "-") ) {
        p_top = stdin;
    } else {
        p_top = fopen(file_name,"r");
        if( p_top == NULL ) {
            CSmallString error;
            error << "unable to open topology file '" << file_name << "' ("
                  << strerror(errno) << ")";
            ES_ERROR(error);
            return(false);
        }
    }

    bool result = LoadFakeTopologyFromPDB(p_top,mangle_names);

    if( result == true ) {
        Name = file_name;
    }

    if( ! ((allow_stdin == true) && (file_name == "-")) ) {
        fclose(p_top);
    }

    return(result);
}

//------------------------------------------------------------------------------

bool CAmberTopology::LoadFakeTopologyFromPDB(FILE* p_fin,bool mangle_names)
{
    Clean();
    FakeTopology = true;

    return(false);
}

//------------------------------------------------------------------------------

bool CAmberTopology::LoadFakeTopologyFromG96(const CSmallString& file_name,bool allow_stdin)
{
    FILE* p_top;

    if( (allow_stdin == true) && (file_name == "-") ) {
        p_top = stdin;
    } else {
        p_top = fopen(file_name,"r");
        if( p_top == NULL ) {
            CSmallString error;
            error << "unable to open topology file '" << file_name << "' ("
                  << strerror(errno) << ")";
            ES_ERROR(error);
            return(false);
        }
    }

    bool result = LoadFakeTopologyFromG96(p_top);

    if( result == true ) {
        Name = file_name;
    }

    if( ! ((allow_stdin == true) && (file_name == "-")) ) {
        fclose(p_top);
    }

    return(result);
}

//------------------------------------------------------------------------------

bool CAmberTopology::LoadFakeTopologyFromG96(FILE* p_fin)
{
    Clean();
    FakeTopology = true;

    list<CFakeResidue> residues;

    CSmallString line;
    while( (feof(p_fin) == 0) && line.ReadLineFromFile(p_fin,false,true) ){
        if( line.GetLength() == 0 ) continue; // empty line
        if( line[0] == '#' ) continue; // comment

        if( line.FindSubString("TITLE") == 0 ){
            line.ReadLineFromFile(p_fin,true,true);
            line.Substitute(' ','_');
            ITITL = line;
        }

    // -----------------
        if( line.FindSubString("POSITION") == 0 ){
            CFakeResidue residue;
            CSmallString lastresid;
            while( (feof(p_fin) == 0) && line.ReadLineFromFile(p_fin,false,true) ){
                if( line.GetLength() == 0 ) continue; // empty line
                if( line[0] == '#' ) continue; // comment
                if( line.FindSubString("END") == 0 ){
                    if( residue.Atoms.size() > 0 ){
                        residues.push_back(residue);
                    }
                    break;
                }

                CSmallString resid = line.GetSubStringFromTo(0,4);
                CSmallString resna = line.GetSubStringFromTo(6,9);
                CSmallString atmna = line.GetSubStringFromTo(12,15);

                if( lastresid == NULL ){
                    residue.Name = resna;
                    lastresid = resid;
                }
                if( lastresid != resid ){
                    residues.push_back(residue);
                    residue.Atoms.clear();
                    residue.Name = resna;
                    lastresid = resid;
                }
                residue.Atoms.push_back(atmna);
            }
            continue;
        }
    // -----------------
        if( line.FindSubString("BOX") == 0 ){
            line.ReadLineFromFile(p_fin,false,true);
            BoxInfo.InitFields(AMBER_BOX_STANDARD);
            BoxInfo.SetBoxBeta(90.0);
            stringstream str(line.GetBuffer());
            CPoint dim;
            str >> dim.x >> dim.y >> dim.z;
            dim *= 10; // nm -> A
            BoxInfo.SetBoxDimmensions(dim);
            continue;
        }
    }

    // determine largest residue
    int maxres = 0;
    list<CFakeResidue>::iterator rit = residues.begin();
    list<CFakeResidue>::iterator rie = residues.end();

    while( rit != rie ){
        int ratoms = (*rit).Atoms.size();
        if( maxres <= ratoms ){
            maxres = ratoms;
        }
        rit++;
    }

    // populate topology

// residues ----
    ResidueList.InitFields(residues.size(),maxres);
    rit = residues.begin();
    int natoms = 0;
    for(int i=0; i < ResidueList.GetNumberOfResidues();i++){
        CAmberResidue* p_res = ResidueList.GetResidue(i);
        CSmallString rname = (*rit).Name;
        int ratoms = (*rit).Atoms.size();
        p_res->SetName(rname);
        p_res->SetFirstAtomIndex(natoms);
        natoms += ratoms;
        rit++;
    }

// atoms ----
    AtomList.InitFields(natoms,0,0);
    rit = residues.begin();
    int i = 0;
    while(rit != rie){
        list<CSmallString>::iterator ait = (*rit).Atoms.begin();
        list<CSmallString>::iterator aie = (*rit).Atoms.end();
        while( ait != aie ){
            CAmberAtom* p_atm = AtomList.GetAtom(i);
            CSmallString aname = (*ait);
            p_atm->SetName(aname);
            i++;
            ait++;
        }
        rit++;
    }

    ResidueList.ReinitAtomResiduePointers(&AtomList);

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberTopology::IsFake(void)
{
    return(FakeTopology);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTopology::Save(const CSmallString& file_name,
                                       EAmberVersion version, bool allow_stdout)
{
    FILE* p_top;

    if( (allow_stdout == true) && (file_name == "-") ) {
        p_top = stdout;
    } else {
        p_top = fopen(file_name,"w");
        if( p_top == NULL ) {
            CSmallString error;
            error << "unable to open topology file '" << file_name << "' ("
                  << strerror(errno) << ")";
            ES_ERROR(error);
            return(false);
        }
    }

    bool result = Save(p_top,version);

    if( ! ((allow_stdout == true) && (file_name == "-")) ) {
        fclose(p_top);
    }

    return(result);
}

//------------------------------------------------------------------------------

bool CAmberTopology::Save(FILE* p_fout,EAmberVersion version)
{
    if( FakeTopology ){
        ES_ERROR("cannot save fake topology");
        return(false);
    }

    bool result;

    if( version == AMBER_VERSION_NONE ) version = Version;

    switch(version) {
    case AMBER_VERSION_6:
        if( Version == AMBER_VERSION_NONE) {
            Version = version;
        }
        result = SaveAmber6(p_fout);
        break;
    case AMBER_VERSION_7:
    default:
        if( Version == AMBER_VERSION_NONE) {
            SetDefaultAmber7Formats();
            Version = version;
        }
        result = SaveAmber7(p_fout);
        break;
    }

    if( result == false ) {
        switch(version) {
        case AMBER_VERSION_6:
            ES_ERROR("unable to save topology file in AMBER 6 format");
            break;
        case AMBER_VERSION_7:
        default:
            ES_ERROR("unable to save topology file in AMBER 7 format");
            break;
        }
    }

    return(result);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTopology::LoadAmber6(FILE* p_top)
{
    if( LoadBasicInfo(p_top,"12I6",AMBER_VERSION_6) == false ) return(false);

    if( AtomList.LoadAtomNames(p_top,"20a4") == false ) return(false);
    if( AtomList.LoadAtomCharges(p_top,"5E16.8") == false ) return(false);
    if( AtomList.LoadAtomMasses(p_top,"5E16.8") == false ) return(false);
    if( AtomList.LoadAtomIACs(p_top,"12I6") == false ) return(false);
    if( AtomList.LoadAtomNUMEXs(p_top,"12I6") == false ) return(false);

    if( NonBondedList.LoadICOs(p_top,"12I6") == false ) return(false);

    if( ResidueList.LoadResidueNames(p_top,"20A4") == false ) return(false);
    if( ResidueList.LoadResidueIPRES(p_top,&AtomList,"12I6") == false ) return(false);

    if( BondList.LoadBondRK(p_top,"5E16.8") == false ) return(false);
    if( BondList.LoadBondREQ(p_top,"5E16.8") == false ) return(false);

    if( AngleList.LoadAngleTK(p_top,"5E16.8") == false ) return(false);
    if( AngleList.LoadAngleTEQ(p_top,"5E16.8") == false ) return(false);

    if( DihedralList.LoadDihedralPK(p_top,"5E16.8") == false ) return(false);
    if( DihedralList.LoadDihedralPN(p_top,"5E16.8") == false ) return(false);
    if( DihedralList.LoadDihedralPHASE(p_top,"5E16.8") == false ) return(false);

    if( NonBondedList.LoadSOLTY(p_top,"5E16.8") == false ) return(false);
    if( NonBondedList.LoadCN1(p_top,"5E16.8") == false ) return(false);
    if( NonBondedList.LoadCN2(p_top,"5E16.8") == false ) return(false);

    if( BondList.LoadBondsWithHydrogens(p_top,"12I6") == false ) return(false);
    if( BondList.LoadBondsWithoutHydrogens(p_top,"12I6") == false ) return(false);

    if( AngleList.LoadAnglesWithHydrogens(p_top,"12I6") == false ) return(false);
    if( AngleList.LoadAnglesWithoutHydrogens(p_top,"12I6") == false ) return(false);

    if( DihedralList.LoadDihedralsWithHydrogens(p_top,"12I6") == false ) return(false);
    if( DihedralList.LoadDihedralsWithoutHydrogens(p_top,"12I6") == false ) return(false);

    if( NonBondedList.LoadNATEX(p_top,"12I6") == false ) return(false);
    if( NonBondedList.LoadASOL(p_top,"5E16.8") == false ) return(false);
    if( NonBondedList.LoadBSOL(p_top,"5E16.8") == false ) return(false);
    if( NonBondedList.LoadHBCUT(p_top,"5E16.8") == false ) return(false);

    if( AtomList.LoadAtomISYMBL(p_top,"20A4") == false ) return(false);
    if( AtomList.LoadAtomITREE(p_top,"20A4") == false ) return(false);
    if( AtomList.LoadAtomJOIN(p_top,"12I6") == false ) return(false);
    if( AtomList.LoadAtomIROTAT(p_top,"12I6") == false ) return(false);

    if( BoxInfo.GetType() != AMBER_BOX_NONE ) { // load box info
        if( BoxInfo.LoadSolventPointers(p_top,"12I6") == false ) return(false);
        if( BoxInfo.LoadNumsOfMolecules(p_top,"12I6") == false ) return(false);
        if( BoxInfo.LoadBoxInfo(p_top,"5E16.8") == false ) return(false);
    }

    if( AtomList.HasPertInfo() == true ) {
        if( BondList.LoadPerturbedBonds(p_top,"12I6") == false ) return(false);
        if( BondList.LoadPerturbedBondTypeIndexes(p_top,"12I6") == false ) return(false);

        if( AngleList.LoadPerturbedAngles(p_top,"12I6") == false ) return(false);
        if( AngleList.LoadPerturbedAngleTypeIndexes(p_top,"12I6") == false ) return(false);

        if( DihedralList.LoadPerturbedDihedrals(p_top,"12I6") == false ) return(false);
        if( DihedralList.LoadPerturbedDihedralTypeIndexes(p_top,"12I6") == false ) return(false);

        if( ResidueList.LoadResiduePertNames(p_top,"20A4") == false ) return(false);

        if( AtomList.LoadPertAtomNames(p_top,"20A4") == false ) return(false);
        if( AtomList.LoadPertAtomISYMBL(p_top,"20A4") == false ) return(false);
        if( AtomList.LoadPertAtomALMPER(p_top,"5E16.8") == false ) return(false);
        if( AtomList.LoadPertAtomPertFlag(p_top,"12I6") == false ) return(false);
        if( AtomList.LoadPertAtomIAC(p_top,"12I6") == false ) return(false);
        if( AtomList.LoadPertAtomCharges(p_top,"5E16.8") == false ) return(false);
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberTopology::SaveAmber6(FILE* p_top)
{
    if( fprintf(p_top,"%-80s\n",(const char*)ITITL) <= 0 ) return(false);

    if( SaveBasicInfo(p_top,"12I6",AMBER_VERSION_6) == false ) return(false);

    if( AtomList.SaveAtomNames(p_top,"20a4") == false ) return(false);
    if( AtomList.SaveAtomCharges(p_top,"5E16.8") == false ) return(false);
    if( AtomList.SaveAtomMasses(p_top,"5E16.8") == false ) return(false);
    if( AtomList.SaveAtomIACs(p_top,"12I6") == false ) return(false);
    if( AtomList.SaveAtomNUMEXs(p_top,"12I6") == false ) return(false);

    if( NonBondedList.SaveICOs(p_top,"12I6") == false ) return(false);

    if( ResidueList.SaveResidueNames(p_top,"20A4") == false ) return(false);
    if( ResidueList.SaveResidueIPRES(p_top,&AtomList,"12I6") == false ) return(false);

    if( BondList.SaveBondRK(p_top,"5E16.8") == false ) return(false);
    if( BondList.SaveBondREQ(p_top,"5E16.8") == false ) return(false);

    if( AngleList.SaveAngleTK(p_top,"5E16.8") == false ) return(false);
    if( AngleList.SaveAngleTEQ(p_top,"5E16.8") == false ) return(false);

    if( DihedralList.SaveDihedralPK(p_top,"5E16.8") == false ) return(false);
    if( DihedralList.SaveDihedralPN(p_top,"5E16.8") == false ) return(false);
    if( DihedralList.SaveDihedralPHASE(p_top,"5E16.8") == false ) return(false);

    if( NonBondedList.SaveSOLTY(p_top,"5E16.8") == false ) return(false);
    if( NonBondedList.SaveCN1(p_top,"5E16.8") == false ) return(false);
    if( NonBondedList.SaveCN2(p_top,"5E16.8") == false ) return(false);

    if( BondList.SaveBondsWithHydrogens(p_top,"12I6") == false ) return(false);
    if( BondList.SaveBondsWithoutHydrogens(p_top,"12I6") == false ) return(false);

    if( AngleList.SaveAnglesWithHydrogens(p_top,"12I6") == false ) return(false);
    if( AngleList.SaveAnglesWithoutHydrogens(p_top,"12I6") == false ) return(false);

    if( DihedralList.SaveDihedralsWithHydrogens(p_top,"12I6") == false ) return(false);
    if( DihedralList.SaveDihedralsWithoutHydrogens(p_top,"12I6") == false ) return(false);

    if( NonBondedList.SaveNATEX(p_top,"12I6") == false ) return(false);
    if( NonBondedList.SaveASOL(p_top,"5E16.8") == false ) return(false);
    if( NonBondedList.SaveBSOL(p_top,"5E16.8") == false ) return(false);
    if( NonBondedList.SaveHBCUT(p_top,"5E16.8") == false ) return(false);

    if( AtomList.SaveAtomISYMBL(p_top,"20A4") == false ) return(false);
    if( AtomList.SaveAtomITREE(p_top,"20A4") == false ) return(false);
    if( AtomList.SaveAtomJOIN(p_top,"12I6") == false ) return(false);
    if( AtomList.SaveAtomIROTAT(p_top,"12I6") == false ) return(false);

    if( BoxInfo.GetType() != AMBER_BOX_NONE ) { // load box info
        if( BoxInfo.SaveSolventPointers(p_top,"12I6") == false ) return(false);
        if( BoxInfo.SaveNumsOfMolecules(p_top,"12I6") == false ) return(false);
        if( BoxInfo.SaveBoxInfo(p_top,"5E16.8") == false ) return(false);
    }

    if( AtomList.HasPertInfo() == true ) {
        if( BondList.SavePerturbedBonds(p_top,"12I6") == false ) return(false);
        if( BondList.SavePerturbedBondTypeIndexes(p_top,"12I6") == false ) return(false);

        if( AngleList.SavePerturbedAngles(p_top,"12I6") == false ) return(false);
        if( AngleList.SavePerturbedAngleTypeIndexes(p_top,"12I6") == false ) return(false);

        if( DihedralList.SavePerturbedDihedrals(p_top,"12I6") == false ) return(false);
        if( DihedralList.SavePerturbedDihedralTypeIndexes(p_top,"12I6") == false ) return(false);

        if( ResidueList.SaveResiduePertNames(p_top,"20A4") == false ) return(false);

        if( AtomList.SavePertAtomNames(p_top,"20A4") == false ) return(false);
        if( AtomList.SavePertAtomISYMBL(p_top,"20A4") == false ) return(false);
        if( AtomList.SavePertAtomALMPER(p_top,"5E16.8") == false ) return(false);
        if( AtomList.SavePertAtomPertFlag(p_top,"12I6") == false ) return(false);
        if( AtomList.SavePertAtomIAC(p_top,"12I6") == false ) return(false);
        if( AtomList.SavePertAtomCharges(p_top,"5E16.8") == false ) return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTopology::LoadAmber7(FILE* p_top)
{
    CFortranIO fortranio(p_top,true);
    char        *p_sname;

    while( (p_sname = fortranio.GetNameOfSection()) != NULL  ) {
        if( strcmp(p_sname,"%FLAG TITLE") == 0 ) {
            fTITLE = fortranio.GetFormatOfSection("%FLAG TITLE");
            if( fTITLE == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG TITLE section");
                return(false);
            }
            // OK - force to read the whole line
            fortranio.SetFormat("1A80");
            if( fortranio.ReadString(ITITL) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG POINTERS") == 0 ) {
            fPOINTERS = fortranio.GetFormatOfSection("%FLAG POINTERS");
            if( fPOINTERS == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG POINTERS section");
                return(false);
            }
            if( LoadBasicInfo(p_top,fPOINTERS,AMBER_VERSION_7) == false ) return(false);
            continue;
        }
        //-----------------------------------
        if( strcmp(p_sname,"%FLAG ATOM_NAME") == 0 ) {
            fATOM_NAME = fortranio.GetFormatOfSection("%FLAG ATOM_NAME");
            if( fATOM_NAME == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG ATOM_NAME section");
                return(false);
            }
            if( AtomList.LoadAtomNames(p_top,fATOM_NAME) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG CHARGE") == 0 ) {
            fCHARGE = fortranio.GetFormatOfSection("%FLAG CHARGE");
            if( fCHARGE == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG CHARGE section");
                return(false);
            }
            if( AtomList.LoadAtomCharges(p_top,fCHARGE) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG ATOMIC_NUMBER") == 0 ) {
            fATOMIC_NUMBER = fortranio.GetFormatOfSection("%FLAG ATOMIC_NUMBER");
            if( fATOMIC_NUMBER == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG ATOMIC_NUMBER section");
                return(false);
            }
            if( AtomList.LoadAtomAtomicNumbers(p_top,fATOMIC_NUMBER) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG MASS") == 0 ) {
            fMASS = fortranio.GetFormatOfSection("%FLAG MASS");
            if( fMASS == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG MASS section");
                return(false);
            }
            if( AtomList.LoadAtomMasses(p_top,fMASS) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG ATOM_TYPE_INDEX") == 0 ) {
            fATOM_TYPE_INDEX = fortranio.GetFormatOfSection("%FLAG ATOM_TYPE_INDEX");
            if( fATOM_TYPE_INDEX == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG ATOM_TYPE_INDEX section");
                return(false);
            }
            if( AtomList.LoadAtomIACs(p_top,fATOM_TYPE_INDEX) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG NUMBER_EXCLUDED_ATOMS") == 0 ) {
            fNUMBER_EXCLUDED_ATOMS = fortranio.GetFormatOfSection("%FLAG NUMBER_EXCLUDED_ATOMS");
            if( fNUMBER_EXCLUDED_ATOMS == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG NUMBER_EXCLUDED_ATOMS section");
                return(false);
            }
            if( AtomList.LoadAtomNUMEXs(p_top,fNUMBER_EXCLUDED_ATOMS) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG NONBONDED_PARM_INDEX") == 0 ) {
            fNONBONDED_PARM_INDEX = fortranio.GetFormatOfSection("%FLAG NONBONDED_PARM_INDEX");
            if( fNONBONDED_PARM_INDEX == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG NONBONDED_PARM_INDEX section");
                return(false);
            }
            if( NonBondedList.LoadICOs(p_top,fNONBONDED_PARM_INDEX) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG RESIDUE_LABEL") == 0 ) {
            fRESIDUE_LABEL = fortranio.GetFormatOfSection("%FLAG RESIDUE_LABEL");
            if( fRESIDUE_LABEL == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG RESIDUE_LABEL section");
                return(false);
            }
            if( ResidueList.LoadResidueNames(p_top,fRESIDUE_LABEL) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG RESIDUE_POINTER") == 0 ) {
            fRESIDUE_POINTER = fortranio.GetFormatOfSection("%FLAG RESIDUE_POINTER");
            if( fRESIDUE_POINTER == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG RESIDUE_POINTER section");
                return(false);
            }
            if( ResidueList.LoadResidueIPRES(p_top,&AtomList,fRESIDUE_POINTER) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG BOND_FORCE_CONSTANT") == 0 ) {
            fBOND_FORCE_CONSTANT = fortranio.GetFormatOfSection("%FLAG BOND_FORCE_CONSTANT");
            if( fBOND_FORCE_CONSTANT == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG BOND_FORCE_CONSTANT section");
                return(false);
            }
            if( BondList.LoadBondRK(p_top,fBOND_FORCE_CONSTANT) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG BOND_EQUIL_VALUE") == 0 ) {
            fBOND_EQUIL_VALUE = fortranio.GetFormatOfSection("%FLAG BOND_EQUIL_VALUE");
            if( fBOND_EQUIL_VALUE == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG BOND_EQUIL_VALUE section");
                return(false);
            }
            if( BondList.LoadBondREQ(p_top,fBOND_EQUIL_VALUE) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG ANGLE_FORCE_CONSTANT") == 0 ) {
            fANGLE_FORCE_CONSTANT = fortranio.GetFormatOfSection("%FLAG ANGLE_FORCE_CONSTANT");
            if( fANGLE_FORCE_CONSTANT == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG ANGLE_FORCE_CONSTANT section");
                return(false);
            }
            if( AngleList.LoadAngleTK(p_top,fANGLE_FORCE_CONSTANT) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG ANGLE_EQUIL_VALUE") == 0 ) {
            fANGLE_EQUIL_VALUE = fortranio.GetFormatOfSection("%FLAG ANGLE_EQUIL_VALUE");
            if( fANGLE_EQUIL_VALUE == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG ANGLE_EQUIL_VALUE section");
                return(false);
            }
            if( AngleList.LoadAngleTEQ(p_top,fANGLE_EQUIL_VALUE) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG DIHEDRAL_FORCE_CONSTANT") == 0 ) {
            fDIHEDRAL_FORCE_CONSTANT = fortranio.GetFormatOfSection("%FLAG DIHEDRAL_FORCE_CONSTANT");
            if( fDIHEDRAL_FORCE_CONSTANT == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG DIHEDRAL_FORCE_CONSTANT section");
                return(false);
            }
            if( DihedralList.LoadDihedralPK(p_top,fDIHEDRAL_FORCE_CONSTANT) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG DIHEDRAL_PERIODICITY") == 0 ) {
            fDIHEDRAL_PERIODICITY = fortranio.GetFormatOfSection("%FLAG DIHEDRAL_PERIODICITY");
            if( fDIHEDRAL_PERIODICITY == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG DIHEDRAL_PERIODICITY section");
                return(false);
            }
            if( DihedralList.LoadDihedralPN(p_top,fDIHEDRAL_PERIODICITY) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG DIHEDRAL_PHASE") == 0 ) {
            fDIHEDRAL_PHASE = fortranio.GetFormatOfSection("%FLAG DIHEDRAL_PHASE");
            if( fDIHEDRAL_PHASE == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG DIHEDRAL_PHASE section");
                return(false);
            }
            if( DihedralList.LoadDihedralPHASE(p_top,fDIHEDRAL_PHASE) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG SCEE_SCALE_FACTOR") == 0 ) {
            fSCEE_SCALE_FACTOR = fortranio.GetFormatOfSection("%FLAG SCEE_SCALE_FACTOR");
            if( fSCEE_SCALE_FACTOR == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG SCEE_SCALE_FACTOR section");
                return(false);
            }
            if( DihedralList.LoadDihedralSCEE(p_top,fSCEE_SCALE_FACTOR) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG SCNB_SCALE_FACTOR") == 0 ) {
            fSCNB_SCALE_FACTOR = fortranio.GetFormatOfSection("%FLAG SCNB_SCALE_FACTOR");
            if( fSCNB_SCALE_FACTOR == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG SCNB_SCALE_FACTOR section");
                return(false);
            }
            if( DihedralList.LoadDihedralSCNB(p_top,fSCNB_SCALE_FACTOR) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG SOLTY") == 0 ) {
            fSOLTY = fortranio.GetFormatOfSection("%FLAG SOLTY");
            if( fSOLTY == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG SOLTY section");
                return(false);
            }
            if( NonBondedList.LoadSOLTY(p_top,fSOLTY) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG LENNARD_JONES_ACOEF") == 0 ) {
            fLENNARD_JONES_ACOEF = fortranio.GetFormatOfSection("%FLAG LENNARD_JONES_ACOEF");
            if( fLENNARD_JONES_ACOEF == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG LENNARD_JONES_ACOEF section");
                return(false);
            }
            if( NonBondedList.LoadCN1(p_top,fLENNARD_JONES_ACOEF) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG LENNARD_JONES_BCOEF") == 0 ) {
            fLENNARD_JONES_BCOEF = fortranio.GetFormatOfSection("%FLAG LENNARD_JONES_BCOEF");
            if( fLENNARD_JONES_BCOEF == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG LENNARD_JONES_BCOEF section");
                return(false);
            }
            if( NonBondedList.LoadCN2(p_top,fLENNARD_JONES_BCOEF) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG BONDS_INC_HYDROGEN") == 0 ) {
            fBONDS_INC_HYDROGEN = fortranio.GetFormatOfSection("%FLAG BONDS_INC_HYDROGEN");
            if( fBONDS_INC_HYDROGEN == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG BONDS_INC_HYDROGEN section");
                return(false);
            }
            if( BondList.LoadBondsWithHydrogens(p_top,fBONDS_INC_HYDROGEN) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG BONDS_WITHOUT_HYDROGEN") == 0 ) {
            fBONDS_WITHOUT_HYDROGEN = fortranio.GetFormatOfSection("%FLAG BONDS_WITHOUT_HYDROGEN");
            if( fBONDS_WITHOUT_HYDROGEN == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG BONDS_WITHOUT_HYDROGEN section");
                return(false);
            }
            if( BondList.LoadBondsWithoutHydrogens(p_top,fBONDS_WITHOUT_HYDROGEN) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG ANGLES_INC_HYDROGEN") == 0 ) {
            fANGLES_INC_HYDROGEN = fortranio.GetFormatOfSection("%FLAG ANGLES_INC_HYDROGEN");
            if( fANGLES_INC_HYDROGEN == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG ANGLES_INC_HYDROGEN section");
                return(false);
            }
            if( AngleList.LoadAnglesWithHydrogens(p_top,fANGLES_INC_HYDROGEN) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG ANGLES_WITHOUT_HYDROGEN") == 0 ) {
            fANGLES_WITHOUT_HYDROGEN = fortranio.GetFormatOfSection("%FLAG ANGLES_WITHOUT_HYDROGEN");
            if( fANGLES_WITHOUT_HYDROGEN == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG ANGLES_WITHOUT_HYDROGEN section");
                return(false);
            }
            if( AngleList.LoadAnglesWithoutHydrogens(p_top,fANGLES_WITHOUT_HYDROGEN) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG DIHEDRALS_INC_HYDROGEN") == 0 ) {
            fDIHEDRALS_INC_HYDROGEN = fortranio.GetFormatOfSection("%FLAG DIHEDRALS_INC_HYDROGEN");
            if( fDIHEDRALS_INC_HYDROGEN == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG DIHEDRALS_INC_HYDROGEN section");
                return(false);
            }
            if( DihedralList.LoadDihedralsWithHydrogens(p_top,fDIHEDRALS_INC_HYDROGEN) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG DIHEDRALS_WITHOUT_HYDROGEN") == 0 ) {
            fDIHEDRALS_WITHOUT_HYDROGEN = fortranio.GetFormatOfSection("%FLAG DIHEDRALS_WITHOUT_HYDROGEN");
            if( fDIHEDRALS_WITHOUT_HYDROGEN == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG DIHEDRALS_WITHOUT_HYDROGEN section");
                return(false);
            }
            if( DihedralList.LoadDihedralsWithoutHydrogens(p_top,fDIHEDRALS_WITHOUT_HYDROGEN) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG EXCLUDED_ATOMS_LIST") == 0 ) {
            fEXCLUDED_ATOMS_LIST = fortranio.GetFormatOfSection("%FLAG EXCLUDED_ATOMS_LIST");
            if( fEXCLUDED_ATOMS_LIST == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG EXCLUDED_ATOMS_LIST section");
                return(false);
            }
            if( NonBondedList.LoadNATEX(p_top,fEXCLUDED_ATOMS_LIST) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG HBOND_ACOEF") == 0 ) {
            fHBOND_ACOEF = fortranio.GetFormatOfSection("%FLAG HBOND_ACOEF");
            if( fHBOND_ACOEF == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG HBOND_ACOEF section");
                return(false);
            }
            if( NonBondedList.LoadASOL(p_top,fHBOND_ACOEF) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG HBOND_BCOEF") == 0 ) {
            fHBOND_BCOEF = fortranio.GetFormatOfSection("%FLAG HBOND_BCOEF");
            if( fHBOND_BCOEF == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG HBOND_BCOEF section");
                return(false);
            }
            if( NonBondedList.LoadBSOL(p_top,fHBOND_BCOEF) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG HBCUT") == 0 ) {
            fHBCUT = fortranio.GetFormatOfSection("%FLAG HBCUT");
            if( fHBCUT == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG HBCUT section");
                return(false);
            }
            if( NonBondedList.LoadHBCUT(p_top,fHBCUT) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG AMBER_ATOM_TYPE") == 0 ) {
            fAMBER_ATOM_TYPE = fortranio.GetFormatOfSection("%FLAG AMBER_ATOM_TYPE");
            if( fAMBER_ATOM_TYPE == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG AMBER_ATOM_TYPE section");
                return(false);
            }
            if( AtomList.LoadAtomISYMBL(p_top,fAMBER_ATOM_TYPE) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG TREE_CHAIN_CLASSIFICATION") == 0 ) {
            fTREE_CHAIN_CLASSIFICATION = fortranio.GetFormatOfSection("%FLAG TREE_CHAIN_CLASSIFICATION");
            if( fTREE_CHAIN_CLASSIFICATION == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG TREE_CHAIN_CLASSIFICATION section");
                return(false);
            }
            if( AtomList.LoadAtomITREE(p_top,fTREE_CHAIN_CLASSIFICATION) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG JOIN_ARRAY") == 0 ) {
            fJOIN_ARRAY = fortranio.GetFormatOfSection("%FLAG JOIN_ARRAY");
            if( fJOIN_ARRAY == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG JOIN_ARRAY section");
                return(false);
            }
            if( AtomList.LoadAtomJOIN(p_top,fJOIN_ARRAY) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG IROTAT") == 0 ) {
            fIROTAT = fortranio.GetFormatOfSection("%FLAG IROTAT");
            if( fIROTAT == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG IROTAT section");
                return(false);
            }
            if( AtomList.LoadAtomIROTAT(p_top,fIROTAT) == false ) return(false);
            continue;
        }

        if( BoxInfo.GetType() != AMBER_BOX_NONE ) { // load box info

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG SOLVENT_POINTERS") == 0 ) {
                fSOLVENT_POINTERS = fortranio.GetFormatOfSection("%FLAG SOLVENT_POINTERS");
                if( fSOLVENT_POINTERS == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG SOLVENT_POINTERS section");
                    return(false);
                }
                if( BoxInfo.LoadSolventPointers(p_top,fSOLVENT_POINTERS) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG ATOMS_PER_MOLECULE") == 0 ) {
                fATOMS_PER_MOLECULE = fortranio.GetFormatOfSection("%FLAG ATOMS_PER_MOLECULE");
                if( fATOMS_PER_MOLECULE == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG ATOMS_PER_MOLECULE section");
                    return(false);
                }
                if( BoxInfo.LoadNumsOfMolecules(p_top,fATOMS_PER_MOLECULE) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG BOX_DIMENSIONS") == 0 ) {
                fBOX_DIMENSIONS = fortranio.GetFormatOfSection("%FLAG BOX_DIMENSIONS");
                if( fBOX_DIMENSIONS == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG BOX_DIMENSIONS section");
                    return(false);
                }
                if( BoxInfo.LoadBoxInfo(p_top,fBOX_DIMENSIONS) == false ) return(false);
                continue;
            }
        }

        //-----------------------------------
        // this is a new section in amber9
        if( strcmp(p_sname,"%FLAG RADIUS_SET") == 0 ) {
            fRADIUS_SET = fortranio.GetFormatOfSection("%FLAG RADIUS_SET");
            if( fRADIUS_SET == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG RADIUS_SET section");
                return(false);
            }
            if( AtomList.LoadAtomRadiusSet(p_top,fRADIUS_SET) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG RADII") == 0 ) {
            fRADII = fortranio.GetFormatOfSection("%FLAG RADII");
            if( fRADII == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG RADII section");
                return(false);
            }
            if( AtomList.LoadAtomRadii(p_top,fRADII) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG SCREEN") == 0 ) {
            fSCREEN = fortranio.GetFormatOfSection("%FLAG SCREEN");
            if( fSCREEN == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG SCREEN section");
                return(false);
            }
            if( AtomList.LoadAtomScreen(p_top,fSCREEN) == false ) return(false);
            continue;
        }


        if( AtomList.HasPertInfo() == true ) {
            //-----------------------------------
            if( strcmp(p_sname,"%FLAG PERT_BOND_ATOMS") == 0 ) {
                fPERT_BOND_ATOMS = fortranio.GetFormatOfSection("%FLAG PERT_BOND_ATOMS");
                if( fPERT_BOND_ATOMS == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG PERT_BOND_ATOMS section");
                    return(false);
                }
                if( BondList.LoadPerturbedBonds(p_top,fPERT_BOND_ATOMS) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG PERT_BOND_PARAMS") == 0 ) {
                fPERT_BOND_PARAMS = fortranio.GetFormatOfSection("%FLAG PERT_BOND_PARAMS");
                if( fPERT_BOND_PARAMS == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG PERT_BOND_PARAMS section");
                    return(false);
                }
                if( BondList.LoadPerturbedBondTypeIndexes(p_top,fPERT_BOND_PARAMS) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG PERT_ANGLE_ATOMS") == 0 ) {
                fPERT_ANGLE_ATOMS = fortranio.GetFormatOfSection("%FLAG PERT_ANGLE_ATOMS");
                if( fPERT_ANGLE_ATOMS == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG PERT_ANGLE_ATOMS section");
                    return(false);
                }
                if( AngleList.LoadPerturbedAngles(p_top,fPERT_ANGLE_ATOMS) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG PERT_ANGLE_PARAMS") == 0 ) {
                fPERT_ANGLE_PARAMS = fortranio.GetFormatOfSection("%FLAG PERT_ANGLE_PARAMS");
                if( fPERT_ANGLE_PARAMS == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG PERT_ANGLE_PARAMS section");
                    return(false);
                }
                if( AngleList.LoadPerturbedAngleTypeIndexes(p_top,fPERT_ANGLE_PARAMS) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG PERT_DIHEDRAL_ATOMS") == 0 ) {
                fPERT_DIHEDRAL_ATOMS = fortranio.GetFormatOfSection("%FLAG PERT_DIHEDRAL_ATOMS");
                if( fPERT_DIHEDRAL_ATOMS == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG PERT_DIHEDRAL_ATOMS section");
                    return(false);
                }
                if( DihedralList.LoadPerturbedDihedrals(p_top,fPERT_DIHEDRAL_ATOMS) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG PERT_DIHEDRAL_PARAMS") == 0 ) {
                fPERT_DIHEDRAL_PARAMS = fortranio.GetFormatOfSection("%FLAG PERT_DIHEDRAL_PARAMS");
                if( fPERT_DIHEDRAL_PARAMS == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG PERT_DIHEDRAL_PARAMS section");
                    return(false);
                }
                if( DihedralList.LoadPerturbedDihedralTypeIndexes(p_top,fPERT_DIHEDRAL_PARAMS) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG PERT_RESIDUE_NAME") == 0 ) {
                fPERT_RESIDUE_NAME = fortranio.GetFormatOfSection("%FLAG PERT_RESIDUE_NAME");
                if( fPERT_RESIDUE_NAME == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG PERT_RESIDUE_NAME section");
                    return(false);
                }
                if( ResidueList.LoadResiduePertNames(p_top,fPERT_RESIDUE_NAME) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG PERT_ATOM_NAME") == 0 ) {
                fPERT_ATOM_NAME = fortranio.GetFormatOfSection("%FLAG PERT_ATOM_NAME");
                if( fPERT_ATOM_NAME == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG PERT_ATOM_NAME section");
                    return(false);
                }
                if( AtomList.LoadPertAtomNames(p_top,fPERT_ATOM_NAME) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG PERT_ATOM_SYMBOL") == 0 ) {
                fPERT_ATOM_SYMBOL = fortranio.GetFormatOfSection("%FLAG PERT_ATOM_SYMBOL");
                if( fPERT_ATOM_SYMBOL == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG PERT_ATOM_SYMBOL section");
                    return(false);
                }
                if( AtomList.LoadPertAtomISYMBL(p_top,fPERT_ATOM_SYMBOL) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG ALMPER") == 0 ) {
                fALMPER = fortranio.GetFormatOfSection("%FLAG ALMPER");
                if( fALMPER == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG ALMPER section");
                    return(false);
                }
                if( AtomList.LoadPertAtomALMPER(p_top,fALMPER) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG IAPER") == 0 ) {
                fIAPER = fortranio.GetFormatOfSection("%FLAG IAPER");
                if( fIAPER == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG IAPER section");
                    return(false);
                }
                if( AtomList.LoadPertAtomPertFlag(p_top,fIAPER) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG  PERT_ATOM_TYPE_INDEX") == 0 ) {
                fPERT_ATOM_TYPE_INDEX = fortranio.GetFormatOfSection("%FLAG  PERT_ATOM_TYPE_INDEX");
                if( fPERT_ATOM_TYPE_INDEX == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG  PERT_ATOM_TYPE_INDEX section");
                    return(false);
                }
                if( AtomList.LoadPertAtomIAC(p_top,fPERT_ATOM_TYPE_INDEX) == false ) return(false);
                continue;
            }

            //-----------------------------------
            if( strcmp(p_sname,"%FLAG PERT_CHARGE") == 0 ) {
                fPERT_CHARGE = fortranio.GetFormatOfSection("%FLAG PERT_CHARGE");
                if( fPERT_CHARGE == NULL ) {
                    ES_ERROR("unable to decode data format of %%FLAG PERT_CHARGE section");
                    return(false);
                }
                if( AtomList.LoadPertAtomCharges(p_top,fPERT_CHARGE) == false ) return(false);
                continue;
            }
        }
        // section was not found
        CSmallString warning;
        warning << "unrecognized section in topology '" << p_sname << ";";
        ES_WARNING(warning);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTopology::SaveAmber7(FILE* p_top)
{
    int outputlen;
    if( (outputlen = fprintf(p_top,"%%VERSION  VERSION_STAMP = V0001.000  DATE = 00/00/00  00:00:00")) <= 0 ) {
        ES_ERROR("unable write version stamp");
        return(false);
    }
    for(int i = outputlen; i < 80; i++) fputc(' ',p_top);
    fputc('\n',p_top);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"TITLE",fTITLE) == false ) return(false);
    CSmallString title = ITITL;
    // trim from left
    int last = title.Verify(" ",-1,-1,true);
    if( last != -1 ){
        title = title.GetSubStringFromTo(0,last);
    }
    if( title == NULL ) title = "default_title";
    title.Substitute(' ','_');
    if( (outputlen = fprintf(p_top,"%s",(const char*)title)) <= 0 ) {
        ES_ERROR("unable write title");
        return(false);
    }
    for(int i = outputlen; i < 80; i++) fputc(' ',p_top);
    fputc('\n',p_top);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"POINTERS",fPOINTERS) == false ) return(false);
    if( SaveBasicInfo(p_top,fPOINTERS,AMBER_VERSION_7) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"ATOM_NAME",fATOM_NAME) == false ) return(false);
    if( AtomList.SaveAtomNames(p_top,fATOM_NAME) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"CHARGE",fCHARGE) == false ) return(false);
    if( AtomList.SaveAtomCharges(p_top,fCHARGE) == false ) return(false);

    //-----------------------------------

    if( AtomList.AtomicNumberLoaded ){
        if( SaveSectionHeader(p_top,"ATOMIC_NUMBER",fATOMIC_NUMBER) == false ) return(false);
        if( AtomList.SaveAtomAtomicNumbers(p_top,fATOMIC_NUMBER) == false ) return(false);
    }

    //-----------------------------------

    if( SaveSectionHeader(p_top,"MASS",fMASS) == false ) return(false);
    if( AtomList.SaveAtomMasses(p_top,fMASS) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"ATOM_TYPE_INDEX",fATOM_TYPE_INDEX) == false ) return(false);
    if( AtomList.SaveAtomIACs(p_top,fATOM_TYPE_INDEX) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"NUMBER_EXCLUDED_ATOMS",fNUMBER_EXCLUDED_ATOMS) == false ) return(false);
    if( AtomList.SaveAtomNUMEXs(p_top,fNUMBER_EXCLUDED_ATOMS) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"NONBONDED_PARM_INDEX",fNONBONDED_PARM_INDEX) == false ) return(false);
    if( NonBondedList.SaveICOs(p_top,fNONBONDED_PARM_INDEX) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"RESIDUE_LABEL",fRESIDUE_LABEL) == false ) return(false);
    if( ResidueList.SaveResidueNames(p_top,fRESIDUE_LABEL) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"RESIDUE_POINTER",fRESIDUE_POINTER) == false ) return(false);
    if( ResidueList.SaveResidueIPRES(p_top,&AtomList,fRESIDUE_POINTER) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"BOND_FORCE_CONSTANT",fBOND_FORCE_CONSTANT) == false ) return(false);
    if( BondList.SaveBondRK(p_top,fBOND_FORCE_CONSTANT) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"BOND_EQUIL_VALUE",fBOND_EQUIL_VALUE) == false ) return(false);
    if( BondList.SaveBondREQ(p_top,fBOND_EQUIL_VALUE) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"ANGLE_FORCE_CONSTANT",fANGLE_FORCE_CONSTANT) == false ) return(false);
    if( AngleList.SaveAngleTK(p_top,fANGLE_FORCE_CONSTANT) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"ANGLE_EQUIL_VALUE",fANGLE_EQUIL_VALUE) == false ) return(false);
    if( AngleList.SaveAngleTEQ(p_top,fANGLE_EQUIL_VALUE) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"DIHEDRAL_FORCE_CONSTANT",fDIHEDRAL_FORCE_CONSTANT) == false ) return(false);
    if( DihedralList.SaveDihedralPK(p_top,fDIHEDRAL_FORCE_CONSTANT) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"DIHEDRAL_PERIODICITY",fDIHEDRAL_PERIODICITY) == false ) return(false);
    if( DihedralList.SaveDihedralPN(p_top,fDIHEDRAL_PERIODICITY) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"DIHEDRAL_PHASE",fDIHEDRAL_PHASE) == false ) return(false);
    if( DihedralList.SaveDihedralPHASE(p_top,fDIHEDRAL_PHASE) == false ) return(false);

    //-----------------------------------

    if( DihedralList.SCEEFactorsLoaded ){
        if( SaveSectionHeader(p_top,"SCEE_SCALE_FACTOR",fSCEE_SCALE_FACTOR) == false ) return(false);
        if( DihedralList.SaveDihedralSCEE(p_top,fSCEE_SCALE_FACTOR) == false ) return(false);
    }

    //-----------------------------------

    if( DihedralList.SCNBFactorsLoaded ){
        if( SaveSectionHeader(p_top,"SCNB_SCALE_FACTOR",fSCNB_SCALE_FACTOR) == false ) return(false);
        if( DihedralList.SaveDihedralSCNB(p_top,fSCNB_SCALE_FACTOR) == false ) return(false);
    }

    //-----------------------------------

    if( SaveSectionHeader(p_top,"SOLTY",fSOLTY) == false ) return(false);
    if( NonBondedList.SaveSOLTY(p_top,fSOLTY) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"LENNARD_JONES_ACOEF",fLENNARD_JONES_ACOEF) == false ) return(false);
    if( NonBondedList.SaveCN1(p_top,fLENNARD_JONES_ACOEF) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"LENNARD_JONES_BCOEF",fLENNARD_JONES_BCOEF) == false ) return(false);
    if( NonBondedList.SaveCN2(p_top,fLENNARD_JONES_BCOEF) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"BONDS_INC_HYDROGEN",fBONDS_INC_HYDROGEN) == false ) return(false);
    if( BondList.SaveBondsWithHydrogens(p_top,fBONDS_INC_HYDROGEN) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"BONDS_WITHOUT_HYDROGEN",fBONDS_WITHOUT_HYDROGEN) == false ) return(false);
    if( BondList.SaveBondsWithoutHydrogens(p_top,fBONDS_WITHOUT_HYDROGEN) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"ANGLES_INC_HYDROGEN",fANGLES_INC_HYDROGEN) == false ) return(false);
    if( AngleList.SaveAnglesWithHydrogens(p_top,fANGLES_INC_HYDROGEN) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"ANGLES_WITHOUT_HYDROGEN",fANGLES_WITHOUT_HYDROGEN) == false ) return(false);
    if( AngleList.SaveAnglesWithoutHydrogens(p_top,fANGLES_WITHOUT_HYDROGEN) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"DIHEDRALS_INC_HYDROGEN",fDIHEDRALS_INC_HYDROGEN) == false ) return(false);
    if( DihedralList.SaveDihedralsWithHydrogens(p_top,fDIHEDRALS_INC_HYDROGEN) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"DIHEDRALS_WITHOUT_HYDROGEN",fDIHEDRALS_WITHOUT_HYDROGEN) == false ) return(false);
    if( DihedralList.SaveDihedralsWithoutHydrogens(p_top,fDIHEDRALS_WITHOUT_HYDROGEN) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"EXCLUDED_ATOMS_LIST",fEXCLUDED_ATOMS_LIST) == false ) return(false);
    if( NonBondedList.SaveNATEX(p_top,fEXCLUDED_ATOMS_LIST) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"HBOND_ACOEF",fHBOND_ACOEF) == false ) return(false);
    if( NonBondedList.SaveASOL(p_top,fHBOND_ACOEF) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"HBOND_BCOEF",fHBOND_BCOEF) == false ) return(false);
    if( NonBondedList.SaveBSOL(p_top,fHBOND_BCOEF) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"HBCUT",fHBCUT) == false ) return(false);
    if( NonBondedList.SaveHBCUT(p_top,fHBCUT) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"AMBER_ATOM_TYPE",fAMBER_ATOM_TYPE) == false ) return(false);
    if( AtomList.SaveAtomISYMBL(p_top,fAMBER_ATOM_TYPE) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"TREE_CHAIN_CLASSIFICATION",fTREE_CHAIN_CLASSIFICATION) == false ) return(false);
    if( AtomList.SaveAtomITREE(p_top,fTREE_CHAIN_CLASSIFICATION) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"JOIN_ARRAY",fJOIN_ARRAY) == false ) return(false);
    if( AtomList.SaveAtomJOIN(p_top,fJOIN_ARRAY) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"IROTAT",fIROTAT) == false ) return(false);
    if( AtomList.SaveAtomIROTAT(p_top,fIROTAT) == false ) return(false);

    //-----------------------------------

    if( BoxInfo.GetType() != AMBER_BOX_NONE ) { // load box info

        if( SaveSectionHeader(p_top,"SOLVENT_POINTERS",fSOLVENT_POINTERS) == false ) return(false);
        if( BoxInfo.SaveSolventPointers(p_top,fSOLVENT_POINTERS) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"ATOMS_PER_MOLECULE",fATOMS_PER_MOLECULE) == false ) return(false);
        if( BoxInfo.SaveNumsOfMolecules(p_top,fATOMS_PER_MOLECULE) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"BOX_DIMENSIONS",fBOX_DIMENSIONS) == false ) return(false);
        if( BoxInfo.SaveBoxInfo(p_top,fBOX_DIMENSIONS) == false ) return(false);
    }

    //-----------------------------------

    if( AtomList.RadiusSet != NULL ) {
        if( SaveSectionHeader(p_top,"RADIUS_SET",fRADIUS_SET) == false ) return(false);
        if( AtomList.SaveAtomRadiusSet(p_top,fRADIUS_SET) == false ) return(false);
    }

    if( SaveSectionHeader(p_top,"RADII",fRADII) == false ) return(false);
    if( AtomList.SaveAtomRadii(p_top,fRADII) == false ) return(false);

    //-----------------------------------

    if( SaveSectionHeader(p_top,"SCREEN",fSCREEN) == false ) return(false);
    if( AtomList.SaveAtomScreen(p_top,fSCREEN) == false ) return(false);

    //-----------------------------------

    if( AtomList.HasPertInfo() == true ) {
        if( SaveSectionHeader(p_top,"PERT_BOND_ATOMS",fPERT_BOND_ATOMS) == false ) return(false);
        if( BondList.SavePerturbedBonds(p_top,fPERT_BOND_ATOMS) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"PERT_BOND_PARAMS",fPERT_BOND_PARAMS) == false ) return(false);
        if( BondList.SavePerturbedBondTypeIndexes(p_top,fPERT_BOND_PARAMS) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"PERT_ANGLE_ATOMS",fPERT_ANGLE_ATOMS) == false ) return(false);
        if( AngleList.SavePerturbedAngles(p_top,fPERT_ANGLE_ATOMS) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"PERT_ANGLE_PARAMS",fPERT_ANGLE_PARAMS) == false ) return(false);
        if( AngleList.SavePerturbedAngleTypeIndexes(p_top,fPERT_ANGLE_PARAMS) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"PERT_DIHEDRAL_ATOMS",fPERT_DIHEDRAL_ATOMS) == false ) return(false);
        if( DihedralList.SavePerturbedDihedrals(p_top,fPERT_DIHEDRAL_ATOMS) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"PERT_DIHEDRAL_PARAMS",fPERT_DIHEDRAL_PARAMS) == false ) return(false);
        if( DihedralList.SavePerturbedDihedralTypeIndexes(p_top,fPERT_DIHEDRAL_PARAMS) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"PERT_RESIDUE_NAME",fPERT_RESIDUE_NAME) == false ) return(false);
        if( ResidueList.SaveResiduePertNames(p_top,fPERT_RESIDUE_NAME) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"PERT_ATOM_NAME",fPERT_ATOM_NAME) == false ) return(false);
        if( AtomList.SavePertAtomNames(p_top,fPERT_ATOM_NAME) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"PERT_ATOM_SYMBOL",fPERT_ATOM_SYMBOL) == false ) return(false);
        if( AtomList.SavePertAtomISYMBL(p_top,fPERT_ATOM_SYMBOL) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"ALMPER",fALMPER) == false ) return(false);
        if( AtomList.SavePertAtomALMPER(p_top,fALMPER) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"IAPER",fIAPER) == false ) return(false);
        if( AtomList.SavePertAtomPertFlag(p_top,fIAPER) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top," PERT_ATOM_TYPE_INDEX",fPERT_ATOM_TYPE_INDEX) == false ) return(false);
        if( AtomList.SavePertAtomIAC(p_top,fPERT_ATOM_TYPE_INDEX) == false ) return(false);

        //-----------------------------------

        if( SaveSectionHeader(p_top,"PERT_CHARGE",fPERT_CHARGE) == false ) return(false);
        if( AtomList.SavePertAtomCharges(p_top,fPERT_CHARGE) == false ) return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTopology::SaveSectionHeader(FILE* p_top,const char* p_section_name,
        const char* p_section_format)
{
    int outputlen;

    if( (outputlen = fprintf(p_top,"%%FLAG %s",p_section_name)) <= 0 ) {
        CSmallString    error;
        error << "unable write header of %%FLAG " << p_section_name << " section";
        ES_ERROR(error);
        return(false);
    }

    for(int i = outputlen; i < 80; i++) fputc(' ',p_top);
    fputc('\n',p_top);

    if( (outputlen = fprintf(p_top,"%%FORMAT(%s)",(char*)p_section_format)) <= 0 ) {
        CSmallString    error;
        error << "unable write format of %%FLAG " << p_section_name << " section";
        ES_ERROR(error);
        return(false);
    }

    for(int i = outputlen; i < 80; i++) fputc(' ',p_top);
    fputc('\n',p_top);

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberTopology::LoadBasicInfo(FILE* p_top,const char* p_format,EAmberVersion version)
{
    CFortranIO fortranio(p_top);
    fortranio.SetFormat(p_format);

    int NATOM;     // total number of atoms
    int NTYPES;    // total number of distinct atom types
    int NBONH;     // number of bonds containing hydrogen
    int MBONA;     // number of bonds not containing hydrogen
    int NTHETH;    // number of angles containing hydrogen
    int MTHETA;    // number of angles not containing hydrogen
    int NPHIH;     // number of dihedrals containing hydrogen
    int MPHIA;     // number of dihedrals not containing hydrogen
    // int NHPARM;    // currently not used  - member item is used
    // int NPARM;     // currently not used  - member item is used
    int NEXT;      // number of excluded atoms
    int NRES;      // number of residues
    int NBONA;     // MBONA + number of constraint bonds
    int NTHETA;    // MTHETA + number of constraint angles
    int NPHIA;     // MPHIA + number of constraint dihedrals
    int NUMBND;    // number of unique bond types
    int NUMANG;    // number of unique angle types
    int NPTRA;     // number of unique dihedral types
    int NATYP;     // number of atom types in parameter file, see SOLTY below
    int NPHB;      // number of distinct 10-12 hydrogen bond pair types
    int IFPERT;    // set to 1 if perturbation info is to be read in
    int NBPER;     // number of bonds to be perturbed
    int NGPER;     // number of angles to be perturbed
    int NDPER;     // number of dihedrals to be perturbed
    int MBPER;     // number of bonds with atoms completely in perturbed group
    int MGPER;     // number of angles with atoms completely in perturbed group
    int MDPER;     // number of dihedrals with atoms completely in perturbed groups
    int IFBOX;     // set to 1 if standard periodic box, 2 when truncated octahedral
    int NMXRS;     // number of atoms in the largest residue
    int IFCAP;     // set to 1 if the CAP option from edit was specified
    int IFPOL=0;     // set to 1 if the polarization option from edit was specified

    if( fortranio.ReadInt(NATOM) == false ) {
        ES_ERROR("unable load NATOM item");
        return(false);
    }

    if( fortranio.ReadInt(NTYPES) == false ) {
        ES_ERROR("unable load NTYPES item");
        return(false);
    }

    if( fortranio.ReadInt(NBONH) == false ) {
        ES_ERROR("unable load NBONH item");
        return(false);
    }

    if( fortranio.ReadInt(MBONA) == false ) {
        ES_ERROR("unable load MBONA item");
        return(false);
    }

    if( fortranio.ReadInt(NTHETH) == false ) {
        ES_ERROR("unable load NTHETH item");
        return(false);
    }

    if( fortranio.ReadInt(MTHETA) == false ) {
        ES_ERROR("unable load MTHETA item");
        return(false);
    }

    if( fortranio.ReadInt(NPHIH) == false ) {
        ES_ERROR("unable load NPHIH item");
        return(false);
    }

    if( fortranio.ReadInt(MPHIA) == false ) {
        ES_ERROR("unable load MPHIA item");
        return(false);
    }

    if( fortranio.ReadInt(NHPARM) == false ) {
        ES_ERROR("unable load NHPARM item");
        return(false);
    }

    if( fortranio.ReadInt(NPARM) == false ) {
        ES_ERROR("unable load NPARM item");
        return(false);
    }

    if( fortranio.ReadInt(NEXT) == false ) {
        ES_ERROR("unable load NEXT item");
        return(false);
    }

    if( fortranio.ReadInt(NRES) == false ) {
        ES_ERROR("unable load atom NRES item");
        return(false);
    }

    if( fortranio.ReadInt(NBONA) == false ) {
        ES_ERROR("unable load NBONA item");
        return(false);
    }

    if( fortranio.ReadInt(NTHETA) == false ) {
        ES_ERROR("unable load atom NTHETA item");
        return(false);
    }

    if( fortranio.ReadInt(NPHIA) == false ) {
        ES_ERROR("unable load NPHIA item");
        return(false);
    }

    if( fortranio.ReadInt(NUMBND) == false ) {
        ES_ERROR("unable load NUMBND item");
        return(false);
    }

    if( fortranio.ReadInt(NUMANG) == false ) {
        ES_ERROR("unable load NUMANG item");
        return(false);
    }

    if( fortranio.ReadInt(NPTRA) == false ) {
        ES_ERROR("unable load NPTRA item");
        return(false);
    }

    if( fortranio.ReadInt(NATYP) == false ) {
        ES_ERROR("unable load NATYP item");
        return(false);
    }

    if( fortranio.ReadInt(NPHB) == false ) {
        ES_ERROR("unable load NPHB item");
        return(false);
    }

    if( fortranio.ReadInt(IFPERT) == false ) {
        ES_ERROR("unable load IFPERT item");
        return(false);
    }

    if( fortranio.ReadInt(NBPER) == false ) {
        ES_ERROR("unable load NBPER item");
        return(false);
    }

    if( fortranio.ReadInt(NGPER) == false ) {
        ES_ERROR("unable load atom NGPER item");
        return(false);
    }

    if( fortranio.ReadInt(NDPER) == false ) {
        ES_ERROR("unable load NDPER item");
        return(false);
    }

    if( fortranio.ReadInt(MBPER) == false ) {
        ES_ERROR("unable load MBPER item");
        return(false);
    }

    if( fortranio.ReadInt(MGPER) == false ) {
        ES_ERROR("unable load MGPER item");
        return(false);
    }

    if( fortranio.ReadInt(MDPER) == false ) {
        ES_ERROR("unable load MDPER item");
        return(false);
    }

    if( fortranio.ReadInt(IFBOX) == false ) {
        ES_ERROR("unable load IFBOX item");
        return(false);
    }

    if( fortranio.ReadInt(NMXRS) == false ) {
        ES_ERROR("unable load NMXRS item");
        return(false);
    }

    if( fortranio.ReadInt(IFCAP) == false ) {
        ES_ERROR("unable load IFCAP item");
        return(false);
    }

    if( version != AMBER_VERSION_6 ) {
        if( fortranio.ReadInt(IFPOL) == false ) {
            ES_ERROR("unable load IFPOL item");
            return(false);
        }
    }

    // init all fields
    AtomList.InitFields(NATOM,IFPERT,IFPOL);
    ResidueList.InitFields(NRES,NMXRS);
    BondList.InitFields(NBONH,MBONA,NUMBND,NBONA,NBPER,MBPER);
    AngleList.InitFields(NTHETH,MTHETA,NUMANG,NTHETA,NGPER,MGPER);
    DihedralList.InitFields(NPHIH,MPHIA,NPTRA,NPHIA,NDPER,MDPER);
    NonBondedList.InitFields(NTYPES,NEXT,NATYP,NPHB);
    BoxInfo.InitFields(IFBOX);
    CapInfo.InitFields(IFCAP);

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberTopology::SaveBasicInfo(FILE* p_top,const char* p_format,EAmberVersion version)
{
    CFortranIO fortranio(p_top);
    fortranio.SetFormat(p_format);

    if( fortranio.WriteInt(AtomList.NATOM) == false ) {
        ES_ERROR("unable save NATOM item");
        return(false);
    }

    if( fortranio.WriteInt(NonBondedList.NTYPES) == false ) {
        ES_ERROR("unable save NTYPES item");
        return(false);
    }

    if( fortranio.WriteInt(BondList.NBONH) == false ) {
        ES_ERROR("unable save NBONH item");
        return(false);
    }

    if( fortranio.WriteInt(BondList.MBONA) == false ) {
        ES_ERROR("unable save MBONA item");
        return(false);
    }

    if( fortranio.WriteInt(AngleList.NTHETH) == false ) {
        ES_ERROR("unable save NTHETH item");
        return(false);
    }

    if( fortranio.WriteInt(AngleList.MTHETA) == false ) {
        ES_ERROR("unable save MTHETA item");
        return(false);
    }

    if( fortranio.WriteInt(DihedralList.NPHIH) == false ) {
        ES_ERROR("unable save NPHIH item");
        return(false);
    }

    if( fortranio.WriteInt(DihedralList.MPHIA) == false ) {
        ES_ERROR("unable save MPHIA item");
        return(false);
    }

    if( fortranio.WriteInt(NHPARM) == false ) {
        ES_ERROR("unable save NHPARM item");
        return(false);
    }

    if( fortranio.WriteInt(NPARM) == false ) {
        ES_ERROR("unable save NPARM item");
        return(false);
    }

    if( fortranio.WriteInt(NonBondedList.NEXT) == false ) {
        ES_ERROR("unable save NEXT item");
        return(false);
    }

    if( fortranio.WriteInt(ResidueList.NRES) == false ) {
        ES_ERROR("unable save NRES item");
        return(false);
    }

    if( fortranio.WriteInt(BondList.NBONA) == false ) {
        ES_ERROR("unable save NBONA item");
        return(false);
    }

    if( fortranio.WriteInt(AngleList.NTHETA) == false ) {
        ES_ERROR("unable save NTHETA item");
        return(false);
    }

    if( fortranio.WriteInt(DihedralList.NPHIA) == false ) {
        ES_ERROR("unable save NPHIA item");
        return(false);
    }

    if( fortranio.WriteInt(BondList.NUMBND) == false ) {
        ES_ERROR("unable save NUMBND item");
        return(false);
    }

    if( fortranio.WriteInt(AngleList.NUMANG) == false ) {
        ES_ERROR("unable save NUMANG item");
        return(false);
    }

    if( fortranio.WriteInt(DihedralList.NPTRA) == false ) {
        ES_ERROR("unable save NPTRA item");
        return(false);
    }

    if( fortranio.WriteInt(NonBondedList.NATYP) == false ) {
        ES_ERROR("unable save NATYP item");
        return(false);
    }

    if( fortranio.WriteInt(NonBondedList.NPHB) == false ) {
        ES_ERROR("unable save NPHB item");
        return(false);
    }

    if( fortranio.WriteInt(AtomList.IFPERT) == false ) {
        ES_ERROR("unable save IFPERT item");
        return(false);
    }

    if( fortranio.WriteInt(BondList.NBPER) == false ) {
        ES_ERROR("unable save NBPER item");
        return(false);
    }

    if( fortranio.WriteInt(AngleList.NGPER) == false ) {
        ES_ERROR("unable save NGPER item");
        return(false);
    }

    if( fortranio.WriteInt(DihedralList.NDPER) == false ) {
        ES_ERROR("unable save NDPER item");
        return(false);
    }

    if( fortranio.WriteInt(BondList.MBPER) == false ) {
        ES_ERROR("unable save MBPER item");
        return(false);
    }

    if( fortranio.WriteInt(AngleList.MGPER) == false ) {
        ES_ERROR("unable save MGPER item");
        return(false);
    }

    if( fortranio.WriteInt(DihedralList.MDPER) == false ) {
        ES_ERROR("unable save MDPER item");
        return(false);
    }

    if( fortranio.WriteInt(BoxInfo.IFBOX) == false ) {
        ES_ERROR("unable save IFBOX item");
        return(false);
    }

    if( fortranio.WriteInt(ResidueList.NMXRS) == false ) {
        ES_ERROR("unable save NMXRS item");
        return(false);
    }

    if( fortranio.WriteInt(CapInfo.IFCAP) == false ) {
        ES_ERROR("unable save IFCAP item");
        return(false);
    }

    if( version != AMBER_VERSION_6 ) {
        if( fortranio.WriteInt(AtomList.IFPOL) == false ) {
            ES_ERROR("unable save IFPOL item");
            return(false);
        }
    }

    fortranio.WriteEndOfSection();

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAmberTopology::SetDefaultAmber7Formats(void)
{
    // amber_7 defaults
    fTITLE="20a4";
    fPOINTERS="10I8";
    fATOM_NAME="20a4";
    fCHARGE="5E16.8";
    fMASS="5E16.8";
    fATOM_TYPE_INDEX="10I8";
    fNUMBER_EXCLUDED_ATOMS="10I8";
    fNONBONDED_PARM_INDEX="10I8";
    fRESIDUE_LABEL="20a4";
    fRESIDUE_POINTER="10I8";
    fBOND_FORCE_CONSTANT="5E16.8";
    fBOND_EQUIL_VALUE="5E16.8";
    fANGLE_FORCE_CONSTANT="5E16.8";
    fANGLE_EQUIL_VALUE="5E16.8";
    fDIHEDRAL_FORCE_CONSTANT="5E16.8";
    fDIHEDRAL_PERIODICITY="5E16.8";
    fDIHEDRAL_PHASE="5E16.8";
    fSOLTY="5E16.8";
    fLENNARD_JONES_ACOEF="5E16.8";
    fLENNARD_JONES_BCOEF="5E16.8";
    fBONDS_INC_HYDROGEN="10I8";
    fBONDS_WITHOUT_HYDROGEN="10I8";
    fANGLES_INC_HYDROGEN="10I8";
    fANGLES_WITHOUT_HYDROGEN="10I8";
    fDIHEDRALS_INC_HYDROGEN="10I8";
    fDIHEDRALS_WITHOUT_HYDROGEN="10I8";
    fEXCLUDED_ATOMS_LIST="10I8";
    fHBOND_ACOEF="5E16.8";
    fHBOND_BCOEF="5E16.8";
    fHBCUT="5E16.8";
    fAMBER_ATOM_TYPE="20a4";
    fTREE_CHAIN_CLASSIFICATION="20a4";
    fJOIN_ARRAY="10I8";
    fIROTAT="10I8";
    fSOLVENT_POINTERS="10I8";
    fATOMS_PER_MOLECULE="10I8";
    fBOX_DIMENSIONS="5E16.8";
    fRADIUS_SET="1A80";
    fRADII="5E16.8";
    fSCREEN="5E16.8";
    fPERT_BOND_ATOMS="";
    fPERT_BOND_PARAMS="";
    fPERT_ANGLE_ATOMS="";
    fPERT_ANGLE_PARAMS="";
    fPERT_DIHEDRAL_ATOMS="";
    fPERT_DIHEDRAL_PARAMS="";
    fPERT_RESIDUE_NAME="";
    fPERT_ATOM_NAME="";
    fPERT_ATOM_SYMBOL="";
    fALMPER="";
    fIAPER="";
    fPERT_ATOM_TYPE_INDEX="";
    fPERT_CHARGE="";
    fSCEE_SCALE_FACTOR = "5E16.8";
    fSCNB_SCALE_FACTOR = "5E16.8";
    fATOMIC_NUMBER = "10I8";
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

