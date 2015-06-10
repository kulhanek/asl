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
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <ErrorSystem.hpp>

#include <AmberTopology.hpp>
#include <AmberMaskRSelection.hpp>
#include <AmberMaskResidues.hpp>
#include <AmberRestart.hpp>

#include "maskparser/AmberMaskParser.hpp"


//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberMaskRSelection::CAmberMaskRSelection(CAmberMaskResidues* p_owner)
{
    Owner = p_owner;
    if( Owner->GetNumberOfTopologyResidues() > 0 ) {
        Residues = new CAmberResidue*[Owner->GetNumberOfTopologyResidues()];
        for(int i=0; i < Owner->GetNumberOfTopologyResidues(); i++) {
            Residues[i] = NULL;
        }
    } else {
        Residues = NULL;
    }
}

//------------------------------------------------------------------------------

CAmberMaskRSelection::~CAmberMaskRSelection(void)
{
    delete[] Residues;
    Residues = NULL;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberMaskRSelection::ExpandAndReduceTree(struct SExpression* p_expr)
{
    if( Residues == NULL ) {
        ES_ERROR("Residues is NULL");
        return(false);
    }

    // expand and reduce tree
    return(ExpandAndReduceTree(this,p_expr));
}

//------------------------------------------------------------------------------

bool CAmberMaskRSelection::ExpandAndReduceTree(
    CAmberMaskRSelection* p_root,
    struct SExpression* p_expr)
{
    if( p_root == NULL ) {
        ES_ERROR("p_root is NULL");
        return(false);
    }

    if( p_expr == NULL ) {
        ES_ERROR("p_expr is NULL");
        return(false);
    }

    // is it selection?
    if( p_expr->Selection != NULL ) {
        if( p_root->Select(p_expr->Selection) == false ) {
            ES_ERROR("unable to select");
            return(false);
        }
        return(true);
    }

    CAmberMaskRSelection* p_left  = NULL;
    CAmberMaskRSelection* p_right = NULL;

    // it is operator - allocate left operand
    switch(p_expr->Operator) {
    case O_NOT:
        // no left operand
        break;
    case O_AND:
    case O_OR:
        try {
            p_left = new CAmberMaskRSelection(p_root->Owner);
        } catch(...) {
            ES_ERROR("unable to allocate left operand");
            return(false);
        }
        if( ExpandAndReduceTree(p_left,p_expr->LeftExpression) == false ) {
            ES_ERROR("cannot expand and reduce left operand");
            delete p_left;
            return(false);
        }
        break;
    case O_ALT:
    case O_AGT:
        ES_ERROR("atom based distance operator is not allowed");
        return(false);
    case O_RLT:
    case O_RGT:
        switch(p_expr->Modificator) {
        case D_ORIGIN:
        case D_CBOX:
            // nothing to do
            break;
        case D_LIST:
        case D_COM:
        case D_PLANE:
            // alocate expression
            try {
                p_left = new CAmberMaskRSelection(p_root->Owner);
            } catch(...) {
                ES_ERROR("unable to allocate left operand");
                return(false);
            }
            // expand and reduce
            if( ExpandAndReduceTree(p_left,p_expr->LeftExpression) == false ) {
                ES_ERROR("cannot expand and reduce left operand");
                delete p_left;
                return(false);
            }
            break;
        default:
            ES_ERROR("not implemented distance modificator");
            return(false);
        }
        break;

    default:
        ES_ERROR("not implemented");
        return(false);
    }

    // it is operator - allocate right operand
    switch(p_expr->Operator) {
    case O_NOT:
    case O_AND:
    case O_OR:
        try {
            p_right = new CAmberMaskRSelection(p_root->Owner);
        } catch(...) {
            ES_ERROR("unable to allocate right operand");
            if( p_left != NULL ) delete p_left;
            return(false);
        }
        break;
    case O_RLT:
    case O_RGT:
        // nothing to do
        break;
    default:
        if( p_left != NULL ) delete p_left;
        ES_ERROR("not implemented");
        return(false);
    }

    // expand and reduce right operand
    if( p_right != NULL ) {
        if( ExpandAndReduceTree(p_right,p_expr->RightExpression) == false ) {
            ES_ERROR("cannot expand and reduce right operand");
            if( p_left != NULL ) delete p_left;
            delete p_right;
            return(false);
        }
    }

    // do arithemetic
    switch(p_expr->Operator) {
    case O_NOT:
        for(int i=0; i < p_root->Owner->GetNumberOfTopologyResidues(); i++) {
            if( p_right->Residues[i] == NULL ) {
                p_root->Residues[i] = p_root->Owner->GetTopology()->ResidueList.GetResidue(i);
            } else {
                p_root->Residues[i] = NULL;
            }
        }
        break;
    case O_AND:
        for(int i=0; i < p_root->Owner->GetNumberOfTopologyResidues(); i++) {
            if( (p_left->Residues[i] != NULL) && (p_right->Residues[i] != NULL) ) {
                p_root->Residues[i] = p_root->Owner->GetTopology()->ResidueList.GetResidue(i);
            } else {
                p_root->Residues[i] = NULL;
            }
        }
        break;
    case O_OR:
        for(int i=0; i < p_root->Owner->GetNumberOfTopologyResidues(); i++) {
            if( (p_left->Residues[i] != NULL) || (p_right->Residues[i] != NULL) ) {
                p_root->Residues[i] = p_root->Owner->GetTopology()->ResidueList.GetResidue(i);
            } else {
                p_root->Residues[i] = NULL;
            }
        }
        break;
    case O_RLT:
    case O_RGT: {
        bool result = false;
        switch(p_expr->Modificator) {
        case D_ORIGIN:
            result = p_root->SelectResidueByDistanceFromOrigin(p_expr->Operator,p_expr->Distance);
            break;
        case D_CBOX:
            result = p_root->SelectResidueByDistanceFromCentreOfBox(p_expr->Operator,p_expr->Distance);
            break;
        case D_LIST:
            result = p_root->SelectResidueByDistanceFromList(p_left,p_expr->Operator,p_expr->Distance);
            break;
        case D_COM:
            result = p_root->SelectResidueByDistanceFromCOM(p_left,p_expr->Operator,p_expr->Distance);
            break;
        case D_PLANE:
            result = p_root->SelectResidueByDistanceFromPlane(p_left,p_expr->Operator,p_expr->Distance);
            break;
        default:
            ES_ERROR("not implemented");
            return(false);
        }
        if( result == false ) return(false);
    }
    break;
    case O_ALT:
    case O_AGT:
        ES_ERROR("atom based distance operator is not allowed");
        break;
    default:
        if( p_left != NULL ) delete p_left;
        if( p_right != NULL ) delete p_right;
        ES_ERROR("not implemented");
        return(false);
    }

    // deallocate operands
    if( p_left != NULL ) delete p_left;
    if( p_right != NULL ) delete p_right;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberMaskRSelection::Select(struct SSelection* p_sel)
{
    if( p_sel == NULL ) {
        ES_ERROR("p_sel is NULL");
        return(false);
    }

    // get selection list
    struct SList* p_list = p_sel->Items;
    if( p_list == NULL ) {
        ES_ERROR("p_list is NULL");
        return(false);
    }

    struct SListItem* p_item = p_list->FirstItem;

    while( p_item != NULL ) {
        if( p_item->Index < 0 ) {
            // no matter of selector - this always means all residues
            for(int i=0; i < Owner->GetNumberOfTopologyResidues(); i++) {
                Residues[i] = Owner->GetTopology()->ResidueList.GetResidue(i);
            }
        }
        if( p_item->Index > 0 ) {
            if( p_item->Length >= 1 ) {
                switch(p_sel->Type) {
                case T_RSELECTOR:
                    SelectResidueByIndex(p_item->Index,p_item->Length);
                    break;
                case T_ASELECTOR:
                    ES_ERROR("atom selector is not permmited");
                    return(false);
                case T_TSELECTOR:
                    ES_ERROR("type selector is not permmited");
                    return(false);
                    return(false);
                default:
                    ES_ERROR("unknown selector");
                    return(false);
                }
            } else {
                ES_ERROR("illegal range");
                return(false);
            }
        }

        if( p_item->Index == 0 ) {
            switch(p_sel->Type) {
            case T_RSELECTOR:
                SelectResidueByName(p_item->Name);
                break;
            case T_ASELECTOR:
                ES_ERROR("atom selector is not permmited");
                return(false);
            case T_TSELECTOR:
                ES_ERROR("atom selector is not permmited");
                return(false);
            default:
                ES_ERROR("unknown selector");
                return(false);
            }
        }

        p_item = p_item->NextItem;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CAmberMaskRSelection::SelectResidueByIndex(int index,int length)
{
    for(int i=0; i < length; i++) {
        int rindex = index + i - 1; // transform to index counted from zero

        // check range - this is mandatory because mask can be general
        if( (rindex >= 0) && (rindex < Owner->GetTopology()->ResidueList.GetNumberOfResidues()) ) {
            Residues[rindex] = Owner->GetTopology()->ResidueList.GetResidue(rindex);
        }
    }
}

//------------------------------------------------------------------------------

void CAmberMaskRSelection::SelectResidueByName(const char* p_name)
{
    // wild card '*' is treated on parser level!
    // see bool Select(struct SSelection* p_sel) method

    // here we need to manage only '=' wild card
    // this wild card can be only at the end of p_name

    int search_len = 4;
    // try to find wild card
    for(int i=3; i >= 0; i--) {
        if( p_name[i] == '=' ) {
            search_len = i;
            break;
        }
    }

    for(int i=0; i < Owner->GetTopology()->ResidueList.GetNumberOfResidues(); i++) {
        CAmberResidue* p_res = Owner->GetTopology()->ResidueList.GetResidue(i);

        if( strncmp(p_res->GetName(),p_name,search_len) == 0 ) {
            // select residue
            Residues[i] = p_res;
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberMaskRSelection::SelectResidueByDistanceFromOrigin(
    SOperator dist_oper,double dist)
{
    if( Owner->GetCoordinates() == NULL ) {
        ES_ERROR("distance operator requires coordinates");
        return(false);
    }

    double dist2 = dist*dist;

    for(int i=0; i < Owner->GetTopology()->ResidueList.GetNumberOfResidues(); i++) {
        CAmberResidue* p_res = Owner->GetTopology()->ResidueList.GetResidue(i);

        bool set = false;

        // check any distance
        for(int j = 0; j < p_res->GetNumberOfAtoms(); j++) {
            int aindex = j + p_res->GetFirstAtomIndex();
            CPoint pos = Owner->GetCoordinates()->GetPosition(aindex);
            double ldist2 = Square(pos);
            switch(dist_oper) {
            case O_RLT:
                if( ldist2 < dist2 ) set = true;
                break;
            case O_RGT:
                if( ldist2 > dist2 ) set = true;
                break;
            case O_ALT:
            case O_AGT:
            default:
                ES_ERROR("incorrect operator");
                return(false);
            }
            if( set == true ) break;
        }

        // select residue
        if( set ) Residues[i] = p_res;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskRSelection::SelectResidueByDistanceFromCentreOfBox(
    SOperator dist_oper,double dist)
{
    if( Owner->GetCoordinates() == NULL ) {
        ES_ERROR("distance operator requires coordinates");
        return(false);
    }

    if( Owner->GetTopology()->BoxInfo.GetType() == AMBER_BOX_NONE ) {
        ES_ERROR("cbox requires box");
        return(false);
    }

    double dist2 = dist*dist;
    CPoint cbox = Owner->GetTopology()->BoxInfo.GetBoxCenter();

    for(int i=0; i < Owner->GetTopology()->ResidueList.GetNumberOfResidues(); i++) {
        CAmberResidue* p_res = Owner->GetTopology()->ResidueList.GetResidue(i);

        bool set = false;

        // check any distance
        for(int j = 0; j < p_res->GetNumberOfAtoms(); j++) {
            int aindex = j + p_res->GetFirstAtomIndex();
            CPoint pos = Owner->GetCoordinates()->GetPosition(aindex);
            double ldist2 = Square(pos-cbox);
            switch(dist_oper) {
            case O_RLT:
                if( ldist2 < dist2 ) set = true;
                break;
            case O_RGT:
                if( ldist2 > dist2 ) set = true;
                break;
            case O_ALT:
            case O_AGT:
            default:
                ES_ERROR("incorrect operator");
                return(false);
            }
            if( set == true ) break;
        }

        // select residue
        if( set ) Residues[i] = p_res;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskRSelection::SelectResidueByDistanceFromList(
    CAmberMaskRSelection* p_left,
    SOperator dist_oper,double dist)
{
    if( Owner->GetCoordinates() == NULL ) {
        ES_ERROR("distance operator requires coordinates");
        return(false);
    }

    double dist2 = dist*dist;

    for(int i=0; i < Owner->GetTopology()->ResidueList.GetNumberOfResidues(); i++) {
        CAmberResidue* p_res = Owner->GetTopology()->ResidueList.GetResidue(i);

        bool set = false;

        for(int j = 0; j < p_res->GetNumberOfAtoms(); j++) {
            int aindex = j + p_res->GetFirstAtomIndex();
            CPoint pos1 = Owner->GetCoordinates()->GetPosition(aindex);

            // check any distance
            for(int k=0; k < Owner->GetTopology()->ResidueList.GetNumberOfResidues(); k++) {
                if( p_left->Residues[k] == NULL ) continue;

                for(int l = 0; l < p_left->Residues[k]->GetNumberOfAtoms(); l++) {
                    int bindex = l + p_left->Residues[k]->GetFirstAtomIndex();

                    CPoint pos2 = Owner->GetCoordinates()->GetPosition(bindex);
                    double ldist2 = Square(pos2-pos1);

                    switch(dist_oper) {
                    case O_RLT:
                        if( ldist2 < dist2 ) {
                            set = true;
                        }
                        break;
                    case O_RGT:
                        if( ldist2 > dist2 ) {
                            set = true;
                        }
                        break;
                    case O_ALT:
                    case O_AGT:
                    default:
                        ES_ERROR("incorrect operator");
                        return(false);
                    }
                    if( set == true ) break;
                }
                if( set == true ) break;
            }
            if( set == true ) break;
        }

        // select residue
        if( set ) Residues[i] = p_res;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskRSelection::SelectResidueByDistanceFromCOM(
    CAmberMaskRSelection* p_left,
    SOperator dist_oper,double dist)
{
    if( Owner->GetCoordinates() == NULL ) {
        ES_ERROR("distance operator requires coordinates");
        return(false);
    }

    // calculate centre of mass
    CPoint     com;
    double     tmass = 0.0;

    for(int i=0; i < Owner->GetTopology()->ResidueList.GetNumberOfResidues(); i++) {
        if( p_left->Residues[i] == NULL ) continue;

        CAmberResidue* p_res = Owner->GetTopology()->ResidueList.GetResidue(i);
        for(int j = 0; j < p_res->GetNumberOfAtoms(); j++) {
            int aindex = j + p_res->GetFirstAtomIndex();
            CAmberAtom* p_atom = Owner->GetTopology()->AtomList.GetAtom(aindex);
            double mass = p_atom->GetMass();
            com += Owner->GetCoordinates()->GetPosition(aindex)*mass;
            tmass += mass;
        }
    }

    if( tmass <= 0.0 ) return(true);
    com /= tmass;

    double dist2 = dist*dist;

    for(int i=0; i < Owner->GetTopology()->ResidueList.GetNumberOfResidues(); i++) {
        CAmberResidue* p_res = Owner->GetTopology()->ResidueList.GetResidue(i);

        bool set = false;

        // check any distance
        for(int j = 0; j < p_res->GetNumberOfAtoms(); j++) {
            int aindex = j + p_res->GetFirstAtomIndex();
            CPoint pos = Owner->GetCoordinates()->GetPosition(aindex);
            double ldist2 = Square(pos-com);
            switch(dist_oper) {
            case O_RLT:
                if( ldist2 < dist2 ) set = true;
                break;
            case O_RGT:
                if( ldist2 > dist2 ) set = true;
                break;
            case O_ALT:
            case O_AGT:
            default:
                ES_ERROR("incorrect operator");
                return(false);
            }
            if( set == true ) break;
        }

        // select residue
        if( set ) Residues[i] = p_res;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskRSelection::SelectResidueByDistanceFromPlane(
    CAmberMaskRSelection* p_left,
    SOperator dist_oper,double dist)
{
    if( Owner->GetCoordinates() == NULL ) {
        ES_ERROR("distance operator requires coordinates");
        return(false);
    }

    ES_ERROR("not implemented");


    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAmberMaskRSelection::GetNumberOfSelectedResidues(void)
{
    if( Residues == NULL ) return(0);

    int num_of_select_residues = 0;

    for(int i=0; i < Owner->GetNumberOfTopologyResidues(); i++) {
        if( Residues[i] != NULL ) num_of_select_residues++;
    }

    return(num_of_select_residues);
}

//------------------------------------------------------------------------------

CAmberResidue* CAmberMaskRSelection::GetSelectedResidue(int index)
{
    if( Residues == NULL ) return(NULL);
    return(Residues[index]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


