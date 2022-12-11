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
#include <vector>
#include <AmberSubTopology.hpp>
#include <AmberMaskAtoms.hpp>
#include <ErrorSystem.hpp>

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberSubTopology::CAmberSubTopology(void)
{
    ReportFile = stdout;

    VerboseMode = false;
    IgnoreErrors = false;
    CopyBox = false;

    AtomMapper = NULL;
    CrossBonds = false;
}

//------------------------------------------------------------------------------

CAmberSubTopology::~CAmberSubTopology(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberSubTopology::InitSubTopology(CAmberMaskAtoms* p_mask,
        bool copy_box,
        bool ignore_errors,
        bool verbose)
{
    // destroy previous data in sub topology
    Clean();

    // do local setup
    VerboseMode = verbose;
    IgnoreErrors = ignore_errors;
    CopyBox = copy_box;

    AtomMapper = NULL;
    CrossBonds = false;

    Mask = p_mask;
    if( Mask != NULL ) OldTopology = Mask->GetTopology();

    bool result = PrepareNewTopology();
    return(result);
}

//------------------------------------------------------------------------------

void CAmberSubTopology::SetReportFile(FILE* p_fout)
{
    ReportFile = p_fout;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberSubTopology::PrepareNewTopology(void)
{
    if( Mask == NULL ) {
        ES_ERROR("mask is NULL");
        return(false);
    }

    if( OldTopology == NULL ) {
        ES_ERROR("old topology is NULL");
        return(false);
    }

    bool result;

    result = PrepareAtoms();

    if( result == true ) result = PrepareResidues();
    if( result == true ) result = PrepareNonbondedList();
    if( result == true ) result = PrepareBonds();
    if( result == true ) result = PrepareAngles();
    if( result == true ) result = PrepareDihedrals();

    if( CopyBox == true ) {
        if( result == true ) result = CopyTopologyBox();
    }

    if( AtomMapper != NULL ) {
        delete[] AtomMapper;
        AtomMapper = NULL;
    }

    CSmallString title;

    if( CrossBonds == true ) {
        title = "cut_topology_with_cross_bonds";
    } else {
        title = "cut_topology";
    }

    title += " [";
    title += LibBuildVersion_ASL;
    title += "]";

    SetTitle(title);

    return(result);
}

//---------------------------------------------------------------------------

bool CAmberSubTopology::PrepareAtoms(void)
{
    int iNATOM;
    int iIFPERT;
    int iIFPOL;
    int iNUMEXTRA;

    if( VerboseMode == false ) {
        fprintf(ReportFile,"Preparing atoms ...\n");
        fprintf(ReportFile,"  Number of atoms : %d\n\n",Mask->GetNumberOfSelectedAtoms());
    }

    AtomMapper = new int[OldTopology->AtomList.GetNumberOfAtoms()];
    if( AtomMapper == NULL ) {
        ES_ERROR("unable to allocate atom mapper");
        return(false);
    }

    iNATOM = Mask->GetNumberOfSelectedAtoms();
    iIFPERT = 0;
    if( OldTopology->AtomList.HasPolInfo() == true ) {
        iIFPOL = 1;
    } else {
        iIFPOL = 0;
    }

    // FIXME
    iNUMEXTRA = 0;

    AtomList.InitFields(iNATOM,iIFPERT,iIFPOL,iNUMEXTRA);

    // set atom data --------------------------------------------------------------
    int j = 0;
    AtomList.AtomicNumberLoaded = OldTopology->AtomList.AtomicNumberLoaded;
    for(int i = 0; i < Mask->GetNumberOfTopologyAtoms(); i++) {
        CAmberAtom* p_oldatom = Mask->GetSelectedAtom(i);
        if( p_oldatom != NULL ) {
            CAmberAtom* p_atom = AtomList.GetAtom(j);
            // copy data
            p_atom->SetMass(p_oldatom->GetMass());
            p_atom->SetName(p_oldatom->GetName());
            p_atom->SetTreeName(p_oldatom->GetTreeName());
            p_atom->SetType(p_oldatom->GetType());
            p_atom->SetRadius(p_oldatom->GetRadius());
            p_atom->SetScreenValue(p_oldatom->GetScreenValue());
            p_atom->SetCharge(p_oldatom->GetCharge());
            p_atom->SetPol(p_oldatom->GetPol());
            p_atom->SetAtomicNumber(p_oldatom->GetAtomicNumber());
            AtomMapper[i] = j;
            j++;
        } else {
            AtomMapper[i] = -1;
        }
    }

    AtomList.SetRadiusSet(OldTopology->AtomList.GetRadiusSet());

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberSubTopology::PrepareResidues(void)
{
    // determine number of residues ---------------------------------------------

    vector<int> residues;

    for(int i = 0; i < Mask->GetNumberOfTopologyAtoms(); i++) {
        CAmberAtom* p_atom = Mask->GetSelectedAtom(i);
        if( p_atom != NULL ) {
            bool new_res = true;
            for(unsigned int j=0; j < residues.size(); j++) {
                if( residues[j] == p_atom->GetResidue()->GetIndex() ) {
                    new_res = false;
                    break;
                }
            }
            if( new_res == true ) residues.push_back(p_atom->GetResidue()->GetIndex());
        }
    }

    // printf info about residues -------------------------------------------------
    if( VerboseMode == false ) {
        fprintf(ReportFile,"Preparing residues ...\n");
        fprintf(ReportFile,"  Number of residues : %d\n\n",(int)residues.size());
    }

    // determine number of atoms in the largest residue ---------------------------
    int maxatom = 0;
    for(unsigned int i=0; i < residues.size(); i++) {
        CAmberResidue* p_oldres = OldTopology->ResidueList.GetResidue(i);
        int real_num_of_atoms = 0;
        for(int j=0; j < p_oldres->GetNumberOfAtoms(); j++) {
            if( Mask->IsAtomSelected(j + p_oldres->GetFirstAtomIndex()) == true ) {
                real_num_of_atoms++;
            }
        }
        if( maxatom < real_num_of_atoms ) maxatom = real_num_of_atoms;
    }

    // allocate fields ------------------------------------------------------------
    ResidueList.InitFields(residues.size(),maxatom);

    // copy data ------------------------------------------------------------------
    int first_atom_index = 0;
    for(unsigned int i=0; i < residues.size(); i++) {
        CAmberResidue* p_res = ResidueList.GetResidue(i);
        CAmberResidue* p_oldres = OldTopology->ResidueList.GetResidue(residues[i]);
        p_res->SetName(p_oldres->GetName());
        p_res->SetFirstAtomIndex(first_atom_index);
        int real_num_of_atoms = 0;
        for(int j=0; j < p_oldres->GetNumberOfAtoms(); j++) {
            if( Mask->IsAtomSelected(j + p_oldres->GetFirstAtomIndex()) == true ) {
                first_atom_index++;
                real_num_of_atoms++;
            }
        }
        // check incomplete residue atom selection
        if( real_num_of_atoms != p_oldres->GetNumberOfAtoms() ) {
            fprintf(ReportFile,"  WARNING: Incomplete residue selection for residue :%03d%s!!\n",p_oldres->GetIndex()+1,p_oldres->GetName());
            fprintf(ReportFile,"           Residue contains %d atoms but only %d is(are) selected!!\n\n",p_oldres->GetNumberOfAtoms(),real_num_of_atoms);
        }
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberSubTopology::PrepareNonbondedList(void)
{
    // determine number of distinct atom types ------------------------------------
    vector<int> iactypes;

    int new_atom_index = 0;
    for(int i = 0; i < Mask->GetNumberOfTopologyAtoms(); i++) {
        CAmberAtom* p_atom = Mask->GetSelectedAtom(i);
        if( p_atom != NULL ) {
            CAmberAtom* p_newatom = AtomList.GetAtom(new_atom_index);

            bool new_type = true;
            for(unsigned int j=0; j < iactypes.size(); j++) {
                if( iactypes[j] == p_atom->GetIAC() ) {
                    new_type = false;
                    p_newatom->SetIAC(j+1);
                    break;
                }
            }
            if( new_type == true ) {
                iactypes.push_back(p_atom->GetIAC());
                p_newatom->SetIAC(iactypes.size());
            }

            new_atom_index++;
        }
    }

    // determine number of excluded atoms -----------------------------------------
    int numex_offset = 0;
    int new_numex_offset = 0;
    new_atom_index = 0;
    for(int i=0; i < OldTopology->AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_oldatom = OldTopology->AtomList.GetAtom(i);
        if( Mask->IsAtomSelected(i) == true ) {
            // determine real number of excluded atoms
            int new_numex = 0;
            for(int j = 0; j < p_oldatom->GetNUMEX(); j++) {
                int nat_index = OldTopology->NonBondedList.GetNATEX(j + numex_offset);
                if( (nat_index >= 0) && (Mask->IsAtomSelected(nat_index) == true) ) {
                    new_numex++;
                }
            }
            if( new_numex == 0 ) {
                new_numex = 1;
            }
            CAmberAtom* p_newatom = AtomList.GetAtom(new_atom_index);
            p_newatom->SetNUMEX(new_numex);
            new_numex_offset += new_numex;
            new_atom_index++;
        }
        numex_offset += p_oldatom->GetNUMEX();
    }

    // determine number of hydrogen bond parameters ------------------------------
    int hbond = 0;
    vector<int> hbondtypes;
    int num_types = iactypes.size();
    for(int i = 0; i < num_types; i++) {
        for(int j = 0; j < i; j++) {
            int index = (iactypes[j]-1)*OldTopology->NonBondedList.GetNumberOfTypes() + (iactypes[i]-1);
            if( OldTopology->NonBondedList.GetICOIndex(index) < 0 ) {
                bool new_type = true;
                for(unsigned int k=0; k < hbondtypes.size(); k++) {
                    if( hbondtypes[k] == OldTopology->NonBondedList.GetICOIndex(index) ) {
                        new_type = false;
                        break;
                    }
                }
                if( new_type == true ) {
                    hbondtypes.push_back(OldTopology->NonBondedList.GetICOIndex(index));
                }
            }
        }
    }
    hbond =  hbondtypes.size();

    // printf info about residues -------------------------------------------------
    if( VerboseMode == false ) {
        fprintf(ReportFile,"Preparing nonbonded list ...\n");
        fprintf(ReportFile,"  Number of distinct types      : %d\n",num_types);
        fprintf(ReportFile,"  Number of excluded atoms      : %d\n",new_numex_offset);
        fprintf(ReportFile,"  Number of hydrogen bond types : %d\n\n",hbond);
    }

    // allocate fields ------------------------------------------------------------
    NonBondedList.InitFields(num_types,new_numex_offset,num_types,hbond);

    // reconstruct CN1 and CN2 arrays ---------------------------------------------
    // CN1 and CN2 arrays are lower trim matrix !
    for(int i = 1; i <= num_types; i++) {
        for(int j = 1; j <= i; j++) {
            int old_param_index;
            int new_param_index;
            if( iactypes[i-1] <= iactypes[j-1] ) {
                old_param_index = iactypes[i-1] + (iactypes[j-1]-1)*iactypes[j-1]/2;
            } else {
                old_param_index = iactypes[j-1] + (iactypes[i-1]-1)*iactypes[i-1]/2;
            }
            new_param_index= j + (i-1)*i/2;
            NonBondedList.SetAParam(OldTopology->NonBondedList.GetAParam(old_param_index),new_param_index);
            NonBondedList.SetBParam(OldTopology->NonBondedList.GetBParam(old_param_index),new_param_index);

            int old_iac = (iactypes[j-1]-1)*OldTopology->NonBondedList.GetNumberOfTypes() + (iactypes[i-1]-1);
            int old_ico = OldTopology->NonBondedList.GetICOIndex(old_iac);
            if( old_ico > 0 ) {
                NonBondedList.SetICOIndex((i-1)*num_types+(j-1),new_param_index);
                NonBondedList.SetICOIndex((j-1)*num_types+(i-1),new_param_index);
            } else {
                // find hydrogen bond
                for(int k = 1; k <= hbond; k++) {
                    if( hbondtypes[k-1] == old_ico ) {
                        NonBondedList.SetICOIndex((i-1)*num_types+(j-1),-k);
                        NonBondedList.SetICOIndex((j-1)*num_types+(i-1),-k);
                        break;
                    }
                }
            }
        }
    }

    // reconstruct ASOL and BSOL arrays -------------------------------------------
    for(int i = 1; i <= hbond; i++) {
        NonBondedList.SetAParam(OldTopology->NonBondedList.GetAParam(hbondtypes[i-1]),-i);
        NonBondedList.SetBParam(OldTopology->NonBondedList.GetBParam(hbondtypes[i-1]),-i);
    }

    // reconstruct NATEX list -----------------------------------------------------
    new_numex_offset = 0;
    numex_offset = 0;
    for(int i=0; i < OldTopology->AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_oldatom = OldTopology->AtomList.GetAtom(i);
        if( Mask->IsAtomSelected(i) == true ) {
            int my_excluded = 0;
            for(int j = 0; j < p_oldatom->GetNUMEX(); j++) {
                int nat_index = OldTopology->NonBondedList.GetNATEX(j + numex_offset);
                if( (nat_index >= 0) && Mask->IsAtomSelected(nat_index) == true ) {
                    NonBondedList.SetNATEX(new_numex_offset,AtomMapper[nat_index]);
                    new_numex_offset++;
                    my_excluded++;
                }
            }
            if( my_excluded == 0 ) {
                NonBondedList.SetNATEX(new_numex_offset,-1);
                new_numex_offset++;
                my_excluded++;
            }
        }
        numex_offset += p_oldatom->GetNUMEX();
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberSubTopology::PrepareBonds(void)
{
    // determine number of bonds with hydrogens -----------------------------------
    int numbondh0 = 0;
    int numbondh1 = 0;

    vector<int> bondwithhydrogens;

    for(int i=0; i < OldTopology->BondList.GetNumberOfBondsWithHydrogen(); i++) {
        //check if bond contains any atom from mask
        bool pass0 = false;   // contain all atoms
        bool pass1 = false;   // contain at least one atom
        CAmberBond* p_bond = OldTopology->BondList.GetBondWithHydrogen(i);
        if( Mask->IsAtomSelected(p_bond->GetIB()) == true ) {
            pass1 = true;
            pass0 = true;
        }
        if( Mask->IsAtomSelected(p_bond->GetJB()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( pass0 == true ) {
            numbondh0++;
            bondwithhydrogens.push_back(i);
        }
        if( pass1 == true ) numbondh1++;
    }

    // determine number of bonds without hydrogens ------------------------------
    int numbond0 = 0;
    int numbond1 = 0;

    vector<int> bondwithouthydrogens;

    for(int i=0; i < OldTopology->BondList.GetNumberOfBondsWithoutHydrogen(); i++) {
        //check if bond contains any atom from mask
        bool pass0 = false;   // contain all atoms
        bool pass1 = false;   // contain at least one atom
        CAmberBond* p_bond = OldTopology->BondList.GetBondWithoutHydrogen(i);
        if( Mask->IsAtomSelected(p_bond->GetIB()) == true ) {
            pass1 = true;
            pass0 = true;
        }
        if( Mask->IsAtomSelected(p_bond->GetJB()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( pass0 == true ) {
            numbond0++;
            bondwithouthydrogens.push_back(i);
        }
        if( pass1 == true ) numbond1++;
    }

    vector<int> bondtypeindexes;

    // determine number of unique bond types --------------------------------------
    for(int i=0; i < OldTopology->BondList.GetNumberOfBonds(); i++) {
        //check if bond contains any atom from mask
        bool pass0 = false;   // contain all atoms
        CAmberBond* p_bond = OldTopology->BondList.GetBond(i);
        if( Mask->IsAtomSelected(p_bond->GetIB()) == true ) {
            pass0 = true;
        }
        if( Mask->IsAtomSelected(p_bond->GetJB()) == true ) {
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( pass0 == true ) {
            // check if it is new type
            bool new_type = true;
            for(unsigned int j=0; j < bondtypeindexes.size(); j++) {
                if( bondtypeindexes[j] == p_bond->GetICB() ) {
                    new_type = false;
                    break;
                }
            }
            if( new_type == true ) bondtypeindexes.push_back(p_bond->GetICB());
        }
    }

    // print info about bonds -----------------------------------------------------

    if( VerboseMode == false ) {
        fprintf(ReportFile,"Preparing bonds ...\n");
        fprintf(ReportFile,"  Number of bonds                  : %d\n",numbond0 + numbondh0);
        fprintf(ReportFile,"  Number of bonds with hydrogens   : %d\n",numbondh0);
        fprintf(ReportFile,"  Number of bonds without hydrogens: %d\n",numbond0);
        fprintf(ReportFile,"  Number of bond types             : %d\n",
                (int)bondtypeindexes.size());
    }
    if( (numbond1 != numbond0) || (numbondh1 != numbondh0) ) {
        CrossBonds = true;
        if( VerboseMode == false ) {
            fprintf(ReportFile,"\n");
        if( IgnoreErrors == false ) {
            fprintf(ReportFile,"  ERROR: Mask selection contains cross bonds!\n");
        } else {
            fprintf(ReportFile,"WARNING: Mask selection contains cross bonds!\n");
            ES_WARNING("mask selection contains cross bonds");
        }
            fprintf(ReportFile,"         Number of cross bonds with hydrogens   : %d\n",numbondh1 - numbondh0);
            fprintf(ReportFile,"         Number of cross bonds without hydrogens: %d\n",numbond1 - numbond0);
        }
        if( IgnoreErrors == false ) {
            ES_ERROR("mask selection contains cross bonds");
            if( VerboseMode == false ) fprintf(ReportFile,"\n");
            return(false);
        }
    }
    if( VerboseMode == false ) {
        fprintf(ReportFile,"\n");
    }

    // allocate fields ------------------------------------------------------------
    BondList.InitFields(numbondh0,numbond0,bondtypeindexes.size(),numbond0,0,0);

    // copy bond types ------------------------------------------------------------
    for(unsigned int i=0; i < bondtypeindexes.size(); i++) {
        CAmberBondType* p_bondtype = BondList.GetBondType(i);
        CAmberBondType* p_oldbondtype = OldTopology->BondList.GetBondType(bondtypeindexes[i]);
        p_bondtype->SetRK(p_oldbondtype->GetRK());
        p_bondtype->SetREQ(p_oldbondtype->GetREQ());
    }

    // copy bonds with hydrogens --------------------------------------------------
    for(unsigned int i=0; i < bondwithhydrogens.size(); i++) {
        CAmberBond* p_bond = BondList.GetBondWithHydrogen(i);
        CAmberBond* p_oldbond = OldTopology->BondList.GetBondWithHydrogen(bondwithhydrogens[i]);
        p_bond->SetIB(AtomMapper[p_oldbond->GetIB()]);
        p_bond->SetJB(AtomMapper[p_oldbond->GetJB()]);
        // find bond type
        unsigned int type = 0;
        for(unsigned int j=0; j < bondtypeindexes.size(); j++) {
            if( bondtypeindexes[j] == p_oldbond->GetICB() ) {
                type = j;
                break;
            }
        }
        p_bond->SetICB(type);
        p_bond->SetPCB(0);
    }

    // copy bonds without hydrogens --------------------------------------------------
    for(unsigned int i=0; i < bondwithouthydrogens.size(); i++) {
        CAmberBond* p_bond = BondList.GetBondWithoutHydrogen(i);
        CAmberBond* p_oldbond = OldTopology->BondList.GetBondWithoutHydrogen(bondwithouthydrogens[i]);
        p_bond->SetIB(AtomMapper[p_oldbond->GetIB()]);
        p_bond->SetJB(AtomMapper[p_oldbond->GetJB()]);
        // find bond type
        unsigned int type = 0;
        for(unsigned int j=0; j < bondtypeindexes.size(); j++) {
            if( bondtypeindexes[j] == p_oldbond->GetICB() ) {
                type = j;
                break;
            }
        }
        p_bond->SetICB(type);
        p_bond->SetPCB(0);
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberSubTopology::PrepareAngles(void)
{
    // determine number of angles with hydrogens --------------------------------
    int numangleh0 = 0;
    int numangleh1 = 0;

    vector<int> anglewithhydrogens;

    for(int i=0; i < OldTopology->AngleList.GetNumberOfAnglesWithHydrogen(); i++) {
        //check if angle contains any atom from mask
        bool pass0 = false;   // contain all atoms
        bool pass1 = false;   // contain at least one atom
        CAmberAngle* p_angle = OldTopology->AngleList.GetAngleWithHydrogen(i);
        if( Mask->IsAtomSelected(p_angle->GetIT()) == true ) {
            pass1 = true;
            pass0 = true;
        }
        if( Mask->IsAtomSelected(p_angle->GetJT()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( Mask->IsAtomSelected(p_angle->GetKT()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( pass0 == true ) {
            numangleh0++;
            anglewithhydrogens.push_back(i);
        }
        if( pass1 == true ) numangleh1++;
    }

    // determine number of angles without hydrogens ------------------------------

    int numangle0 = 0;
    int numangle1 = 0;

    vector<int> anglewithouthydrogens;

    for(int i=0; i < OldTopology->AngleList.GetNumberOfAnglesWithoutHydrogen(); i++) {
        //check if angle contains any atom from mask
        bool pass0 = false;   // contain all atoms
        bool pass1 = false;   // contain at least one atom
        CAmberAngle* p_angle = OldTopology->AngleList.GetAngleWithoutHydrogen(i);
        if( Mask->IsAtomSelected(p_angle->GetIT()) == true ) {
            pass1 = true;
            pass0 = true;
        }
        if( Mask->IsAtomSelected(p_angle->GetJT()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( Mask->IsAtomSelected(p_angle->GetKT()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( pass0 == true ) {
            numangle0++;
            anglewithouthydrogens.push_back(i);
        }
        if( pass1 == true ) numangle1++;
    }

    vector<int> angletypeindexes;

    // determine number of unique angle types -------------------------------------

    for(int i=0; i < OldTopology->AngleList.GetNumberOfAngles(); i++) {
        //check if angle contains any atom from mask
        bool pass0 = false;   // contain all atoms
        CAmberAngle* p_angle = OldTopology->AngleList.GetAngle(i);
        if( Mask->IsAtomSelected(p_angle->GetIT()) == true ) {
            pass0 = true;
        }
        if( Mask->IsAtomSelected(p_angle->GetJT()) == true ) {
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( Mask->IsAtomSelected(p_angle->GetKT()) == true ) {
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( pass0 == true ) {
            // check if it is new type
            bool new_type = true;
            for(unsigned int j=0; j < angletypeindexes.size(); j++) {
                if( angletypeindexes[j] == p_angle->GetICT() ) {
                    new_type = false;
                    break;
                }
            }
            if( new_type == true ) angletypeindexes.push_back(p_angle->GetICT());
        }
    }

    // print info about angles -----------------------------------------------------

    if( VerboseMode == false ) {
        fprintf(ReportFile,"Preparing angles ...\n");
        fprintf(ReportFile,"  Number of angles                  : %d\n",numangle0 + numangleh0);
        fprintf(ReportFile,"  Number of angles with hydrogens   : %d\n",numangleh0);
        fprintf(ReportFile,"  Number of angles without hydrogens: %d\n",numangle0);
        fprintf(ReportFile,"  Number of angle types             : %d\n",
                (int)angletypeindexes.size());
    }
    if( (numangle1 != numangle0) || (numangleh1 != numangleh0) ) {
        if( VerboseMode == false ) {
            fprintf(ReportFile,"\n");
        if( IgnoreErrors == false ) {
            fprintf(ReportFile,"  ERROR: Mask selection contains cross angles!\n");
        } else {
            fprintf(ReportFile,"WARNING: Mask selection contains cross angles!\n");
            ES_WARNING("mask selection contains cross angles");
        }
            fprintf(ReportFile,"         Number of cross angles with hydrogens   : %d\n",numangleh1 - numangleh0);
            fprintf(ReportFile,"         Number of cross angles without hydrogens: %d\n",numangle1 - numangle0);
        }
        if( IgnoreErrors == false ) {
            ES_ERROR("mask selection contains cross angles");
            if( VerboseMode == false ) fprintf(ReportFile,"\n");
            return(false);
        }
    }
    if( VerboseMode == false ) {
        fprintf(ReportFile,"\n");
    }

    // allocate fields ------------------------------------------------------------
    AngleList.InitFields(numangleh0,numangle0,angletypeindexes.size(),numangle0,0,0);

    // copy angle types ------------------------------------------------------------
    for(unsigned int i=0; i < angletypeindexes.size(); i++) {
        CAmberAngleType* p_angletype = AngleList.GetAngleType(i);
        CAmberAngleType* p_oldangletype = OldTopology->AngleList.GetAngleType(angletypeindexes[i]);
        p_angletype->SetTK(p_oldangletype->GetTK());
        p_angletype->SetTEQ(p_oldangletype->GetTEQ());
    }

    // copy angles with hydrogens --------------------------------------------------
    for(unsigned int i=0; i < anglewithhydrogens.size(); i++) {
        CAmberAngle* p_angle = AngleList.GetAngleWithHydrogen(i);
        CAmberAngle* p_oldangle = OldTopology->AngleList.GetAngleWithHydrogen(anglewithhydrogens[i]);
        p_angle->SetIT(AtomMapper[p_oldangle->GetIT()]);
        p_angle->SetJT(AtomMapper[p_oldangle->GetJT()]);
        p_angle->SetKT(AtomMapper[p_oldangle->GetKT()]);
        // find bond type
        unsigned int type = 0;
        for(unsigned int j=0; j < angletypeindexes.size(); j++) {
            if( angletypeindexes[j] == p_oldangle->GetICT() ) {
                type = j;
                break;
            }
        }
        p_angle->SetICT(type);
        p_angle->SetPCT(0);
    }

    // copy angles without hydrogens --------------------------------------------------
    for(unsigned int i=0; i < anglewithouthydrogens.size(); i++) {
        CAmberAngle* p_angle = AngleList.GetAngleWithoutHydrogen(i);
        CAmberAngle* p_oldangle = OldTopology->AngleList.GetAngleWithoutHydrogen(anglewithouthydrogens[i]);
        p_angle->SetIT(AtomMapper[p_oldangle->GetIT()]);
        p_angle->SetJT(AtomMapper[p_oldangle->GetJT()]);
        p_angle->SetKT(AtomMapper[p_oldangle->GetKT()]);
        // find bond type
        unsigned int type = 0;
        for(unsigned int j=0; j < angletypeindexes.size(); j++) {
            if( angletypeindexes[j] == p_oldangle->GetICT() ) {
                type = j;
                break;
            }
        }
        p_angle->SetICT(type);
        p_angle->SetPCT(0);
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberSubTopology::PrepareDihedrals(void)
{
    // determine number of dihedrals with hydrogens --------------------------------
    int numdihedralh0 = 0;
    int numdihedralh1 = 0;

    vector<int> dihedralwithhydrogens;

    for(int i=0; i < OldTopology->DihedralList.GetNumberOfDihedralsWithHydrogen(); i++) {
        //check if dihedral contains any atom from mask
        bool pass0 = false;   // contain all atoms
        bool pass1 = false;   // contain at least one atom
        CAmberDihedral* p_dihedral = OldTopology->DihedralList.GetDihedralWithHydrogen(i);
        if( Mask->IsAtomSelected(p_dihedral->GetIP()) == true ) {
            pass1 = true;
            pass0 = true;
        }
        if( Mask->IsAtomSelected(p_dihedral->GetJP()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( Mask->IsAtomSelected(p_dihedral->GetKP()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( Mask->IsAtomSelected(p_dihedral->GetLP()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( pass0 == true ) {
            numdihedralh0++;
            dihedralwithhydrogens.push_back(i);
        }
        if( pass1 == true ) numdihedralh1++;
    }

    // determine number of dihedrals without hydrogens ------------------------------
    int numdihedral0 = 0;
    int numdihedral1 = 0;

    vector<int> dihedralwithouthydrogens;

    for(int i=0; i < OldTopology->DihedralList.GetNumberOfDihedralsWithoutHydrogen(); i++) {
        //check if dihedral contains any atom from mask
        bool pass0 = false;   // contain all atoms
        bool pass1 = false;   // contain at least one atom
        CAmberDihedral* p_dihedral = OldTopology->DihedralList.GetDihedralWithoutHydrogen(i);
        if( Mask->IsAtomSelected(p_dihedral->GetIP()) == true ) {
            pass1 = true;
            pass0 = true;
        }
        if( Mask->IsAtomSelected(p_dihedral->GetJP()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( Mask->IsAtomSelected(p_dihedral->GetKP()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( Mask->IsAtomSelected(p_dihedral->GetLP()) == true ) {
            pass1 = true;
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( pass0 == true ) {
            numdihedral0++;
            dihedralwithouthydrogens.push_back(i);
        }
        if( pass1 == true ) numdihedral1++;
    }

    vector<int> dihedraltypeindexes;

    // determine number of unique dihedral types --------------------------------------
    for(int i=0; i < OldTopology->DihedralList.GetNumberOfDihedrals(); i++) {
        //check if dihedral contains any atom from mask
        bool pass0 = false;   // contain all atoms
        CAmberDihedral* p_dihedral = OldTopology->DihedralList.GetDihedral(i);
        if( Mask->IsAtomSelected(p_dihedral->GetIP()) == true ) {
            pass0 = true;
        }
        if( Mask->IsAtomSelected(p_dihedral->GetJP()) == true ) {
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( Mask->IsAtomSelected(p_dihedral->GetKP()) == true ) {
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( Mask->IsAtomSelected(p_dihedral->GetLP()) == true ) {
            pass0 &= true;
        } else {
            pass0 &= false;
        }
        if( pass0 == true ) {
            // check if it is new type
            bool new_type = true;
            for(unsigned int j=0; j < dihedraltypeindexes.size(); j++) {
                if( dihedraltypeindexes[j] == p_dihedral->GetICP() ) {
                    new_type = false;
                    break;
                }
            }
            if( new_type == true ) dihedraltypeindexes.push_back(p_dihedral->GetICP());
        }
    }

    // print info about dihedrals -----------------------------------------------------

    if( VerboseMode == false ) {
        fprintf(ReportFile,"Preparing dihedrals ...\n");
        fprintf(ReportFile,"  Number of dihedrals                  : %d\n",numdihedral0 + numdihedralh0);
        fprintf(ReportFile,"  Number of dihedrals with hydrogens   : %d\n",numdihedralh0);
        fprintf(ReportFile,"  Number of dihedrals without hydrogens: %d\n",numdihedral0);
        fprintf(ReportFile,"  Number of dihedral types             : %d\n",
                (int)dihedraltypeindexes.size());
    }
    if( (numdihedral1 != numdihedral0) || (numdihedralh1 != numdihedralh0) ) {
        if( VerboseMode == false ) {
            fprintf(ReportFile,"\n");
        if( IgnoreErrors == false ) {
            fprintf(ReportFile,"  ERROR: Mask selection contains cross dihedrals!\n");
        } else {
            fprintf(ReportFile,"WARNING: Mask selection contains cross dihedrals!\n");
            ES_WARNING("mask selection contains cross dihedrals");
        }

            fprintf(ReportFile,"         Number of cross dihedrals with hydrogens   : %d\n",numdihedralh1 - numdihedralh0);
            fprintf(ReportFile,"         Number of cross dihedrals without hydrogens: %d\n",numdihedral1 - numdihedral0);
        }
        if( IgnoreErrors == false ) {
            ES_ERROR("mask selection contains cross dihedrals");
            if( VerboseMode == false ) fprintf(ReportFile,"\n");
            return(false);
        }
    }
    if( VerboseMode == false ) {
        fprintf(ReportFile,"\n");
    }

    // allocate fields ------------------------------------------------------------
    DihedralList.InitFields(numdihedralh0,numdihedral0,dihedraltypeindexes.size(),numdihedral0,0,0);

    // copy dihedral types ------------------------------------------------------------
    DihedralList.SCEEFactorsLoaded = OldTopology->DihedralList.SCEEFactorsLoaded;
    DihedralList.SCNBFactorsLoaded = OldTopology->DihedralList.SCNBFactorsLoaded;
    for(unsigned int i=0; i < dihedraltypeindexes.size(); i++) {
        CAmberDihedralType* p_dihedraltype = DihedralList.GetDihedralType(i);
        CAmberDihedralType* p_olddihedraltype = OldTopology->DihedralList.GetDihedralType(dihedraltypeindexes[i]);
        p_dihedraltype->SetPK(p_olddihedraltype->GetPK());
        p_dihedraltype->SetPN(p_olddihedraltype->GetPN());
        p_dihedraltype->SetPHASE(p_olddihedraltype->GetPHASE());
        p_dihedraltype->SetSCEE(p_olddihedraltype->GetSCEE());
        p_dihedraltype->SetSCNB(p_olddihedraltype->GetSCNB());
    }

    int pokus = 0;

    // copy angles with hydrogens --------------------------------------------------
    for(unsigned int i=0; i < dihedralwithhydrogens.size(); i++) {
        CAmberDihedral* p_dihedral = DihedralList.GetDihedralWithHydrogen(i);
        CAmberDihedral* p_olddihedral= OldTopology->DihedralList.GetDihedralWithHydrogen(dihedralwithhydrogens[i]);

        // WARNING !!!!! : KP and LP value must be non-zero values !!!!!
        //                 becase sign(KP) and sign(LP) determine type of dihedral (sign for zero number is not defined)

        if( (AtomMapper[p_olddihedral->GetKP()] == 0) || (AtomMapper[p_olddihedral->GetLP()] == 0) ) {
            p_dihedral->SetIP(AtomMapper[p_olddihedral->GetLP()]);
            p_dihedral->SetJP(AtomMapper[p_olddihedral->GetKP()]);
            p_dihedral->SetKP(AtomMapper[p_olddihedral->GetJP()]);
            p_dihedral->SetLP(AtomMapper[p_olddihedral->GetIP()]);
            p_dihedral->SetType(p_olddihedral->GetType());
        } else {
            p_dihedral->SetIP(AtomMapper[p_olddihedral->GetIP()]);
            p_dihedral->SetJP(AtomMapper[p_olddihedral->GetJP()]);
            p_dihedral->SetKP(AtomMapper[p_olddihedral->GetKP()]);
            p_dihedral->SetLP(AtomMapper[p_olddihedral->GetLP()]);
            p_dihedral->SetType(p_olddihedral->GetType());
        }
        if( p_olddihedral->GetType() == 0 ) pokus++;

        // find bond type
        unsigned int type = 0;
        for(unsigned int j=0; j < dihedraltypeindexes.size(); j++) {
            if( dihedraltypeindexes[j] == p_olddihedral->GetICP() ) {
                type = j;
                break;
            }
        }
        p_dihedral->SetICP(type);
        p_dihedral->SetPCP(0);
    }

    // copy angles without hydrogens --------------------------------------------------
    for(unsigned int i=0; i < dihedralwithouthydrogens.size(); i++) {
        CAmberDihedral* p_dihedral = DihedralList.GetDihedralWithoutHydrogen(i);
        CAmberDihedral* p_olddihedral = OldTopology->DihedralList.GetDihedralWithoutHydrogen(dihedralwithouthydrogens[i]);

        // WARNING !!!!! : KP and LP value must be non-zero values !!!!!
        //                 becase sign(KP) and sign(LP) determine type of dihedral (sign for zero number is not defined)

        if( ((AtomMapper[p_olddihedral->GetKP()]) == 0) || ((AtomMapper[p_olddihedral->GetKP()]) == 0) ) {
            p_dihedral->SetIP(AtomMapper[p_olddihedral->GetLP()]);
            p_dihedral->SetJP(AtomMapper[p_olddihedral->GetKP()]);
            p_dihedral->SetKP(AtomMapper[p_olddihedral->GetJP()]);
            p_dihedral->SetLP(AtomMapper[p_olddihedral->GetIP()]);
            p_dihedral->SetType(p_olddihedral->GetType());
        } else {
            p_dihedral->SetIP(AtomMapper[p_olddihedral->GetIP()]);
            p_dihedral->SetJP(AtomMapper[p_olddihedral->GetJP()]);
            p_dihedral->SetKP(AtomMapper[p_olddihedral->GetKP()]);
            p_dihedral->SetLP(AtomMapper[p_olddihedral->GetLP()]);
            p_dihedral->SetType(p_olddihedral->GetType());
        }

        // find bond type
        unsigned int type = 0;
        for(unsigned int j=0; j < dihedraltypeindexes.size(); j++) {
            if( dihedraltypeindexes[j] == p_olddihedral->GetICP() ) {
                type = j;
                break;
            }
        }
        p_dihedral->SetICP(type);
        p_dihedral->SetPCP(0);
    }

    return(true);
}

//---------------------------------------------------------------------------

bool CAmberSubTopology::CopyTopologyBox(void)
{
    if( VerboseMode == false ) {
        fprintf(ReportFile,"Preparing box ...\n");
    }

    if( OldTopology->BoxInfo.GetType() == AMBER_BOX_NONE ) {
        if( VerboseMode == false ) {
            fprintf(ReportFile," WARNING: Input topology does not contain box information!\n");
            fprintf(ReportFile,"\n");
        }
        return(true);
    }

    if( VerboseMode == false ) {
        fprintf(ReportFile,"  Box type                  :");
        switch( OldTopology->BoxInfo.GetType() ) {
        case AMBER_BOX_STANDARD:
            fprintf(ReportFile," rectangular\n");
            break;
        case AMBER_BOX_OCTAHEDRAL:
            fprintf(ReportFile," truncated octahedral\n");
            break;
        case AMBER_BOX_NONE:
            // this cannot happen
            return(false);
        }
    }

    if( CrossBonds == true ) {
        if( IgnoreErrors == false ){
            ES_ERROR("selection contains cross bonds - box cannot be copied");
            return(false);
        }
    }

    // copy basic data
    BoxInfo.SetType(OldTopology->BoxInfo.GetType());

    BoxInfo.SetBoxDimmensions(OldTopology->BoxInfo.GetBoxDimmensions());

    BoxInfo.SetBoxBeta(OldTopology->BoxInfo.GetBoxBeta());

    // prepare molecules
    ResidueList.ReinitAtomResiduePointers(&AtomList);
    InitMoleculeIndexes();

    // calculate number of molecules
    int mol_count = 0;
    int prev_index = -1;
    for(int i=0; i < AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = AtomList.GetAtom(i);
        if( prev_index == -1 ) {
            prev_index = p_atom->GetMoleculeIndex();
            mol_count++;
        } else {
            if( prev_index != p_atom->GetMoleculeIndex() ) {
                prev_index = p_atom->GetMoleculeIndex();
                mol_count++;
            }
        }
    }

    if( VerboseMode == false )fprintf(ReportFile,"  Number of molecules: %d\n\n",mol_count);
    BoxInfo.SetNumberOfMolecules(mol_count);

    // set atom numbers in each molecule
    int mol_index = -1;
    prev_index = -1;
    int mol_natoms = 0;
    for(int i=0; i < AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = AtomList.GetAtom(i);
        if( prev_index == -1 ) {
            prev_index = p_atom->GetMoleculeIndex();
            mol_index++;
            mol_natoms++;
        } else {
            if( prev_index != p_atom->GetMoleculeIndex() ) {
                prev_index = p_atom->GetMoleculeIndex();
                BoxInfo.SetNumberOfAtomsInMolecule(mol_index,mol_natoms);
                mol_index++;
                mol_natoms = 1;
            } else {
                mol_natoms++;
            }
        }
    }
    if( mol_index != -1 ) {
        BoxInfo.SetNumberOfAtomsInMolecule(mol_index,mol_natoms);
    }

    // determine first solvent molecule -------------
    // go over list of old solute molecules
    int atom_idx = 0;
    int new_num_of_solute_mols = 0;
    for(int i=0; i < OldTopology->BoxInfo.GetNumberOfSoluteMolecules(); i++) {
        bool keep = false;
        for(int j=atom_idx; j < atom_idx+OldTopology->BoxInfo.GetNumberOfAtomsInMolecule(i); j++) {
            keep |= Mask->IsAtomSelected(j);
        }
        atom_idx += OldTopology->BoxInfo.GetNumberOfAtomsInMolecule(i);
        if( keep ) new_num_of_solute_mols++;
    }

    BoxInfo.SetFirstSolventMolecule(new_num_of_solute_mols);

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
