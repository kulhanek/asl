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
#include <AmberRestart.hpp>
#include <AmberMaskASelection.hpp>
#include <AmberMaskAtoms.hpp>

#include "maskparser/AmberMaskParser.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAmberMaskASelection::CAmberMaskASelection(CAmberMaskAtoms* p_owner)
{
    Owner = p_owner;
    if( Owner->GetNumberOfTopologyAtoms() > 0 ) {
        Atoms = new CAmberAtom*[Owner->GetNumberOfTopologyAtoms()];
        for(int i=0; i < Owner->GetNumberOfTopologyAtoms(); i++) Atoms[i] = NULL;
    } else {
        Atoms = NULL;
    }
}

//------------------------------------------------------------------------------

CAmberMaskASelection::~CAmberMaskASelection(void)
{
    delete[] Atoms;
    Atoms = NULL;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberMaskASelection::ExpandAndReduceTree(struct SExpression* p_expr)
{
    if( Atoms == NULL ) {
        ES_ERROR("Atoms is NULL");
        return(false);
    }

    // expand and reduce tree
    return(ExpandAndReduceTree(this,p_expr));
}

//------------------------------------------------------------------------------

bool CAmberMaskASelection::ExpandAndReduceTree(
    CAmberMaskASelection* p_root,
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

    CAmberMaskASelection* p_left  = NULL;
    CAmberMaskASelection* p_right = NULL;

    // it is operator - allocate left operand
    switch(p_expr->Operator) {
    case O_NOT:
        // no left operand
        break;
    case O_AND:
    case O_OR:
        // alocate expression
        try {
            p_left = new CAmberMaskASelection(p_root->Owner);
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
    case O_RLT:
    case O_RGT:
    case O_ALT:
    case O_AGT:
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
                p_left = new CAmberMaskASelection(p_root->Owner);
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
        ES_ERROR("not implemented operator");
        return(false);
    }

    // it is operator - allocate right operand
    switch(p_expr->Operator) {
    case O_NOT:
    case O_AND:
    case O_OR:
        try {
            p_right = new CAmberMaskASelection(p_root->Owner);
        } catch(...) {
            ES_ERROR("unable to allocate right operand");
            if( p_left != NULL ) delete p_left;
            return(false);
        }
        break;
    case O_RLT:
    case O_RGT:
    case O_ALT:
    case O_AGT:
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
        for(int i=0; i < p_root->Owner->GetNumberOfTopologyAtoms(); i++) {
            if( p_right->Atoms[i] == NULL ) {
                p_root->Atoms[i] = p_root->Owner->GetTopology()->AtomList.GetAtom(i);
            } else {
                p_root->Atoms[i] = NULL;
            }
        }
        break;
    case O_AND:
        for(int i=0; i < p_root->Owner->GetNumberOfTopologyAtoms(); i++) {
            if( (p_left->Atoms[i] != NULL) && (p_right->Atoms[i] != NULL) ) {
                p_root->Atoms[i] = p_root->Owner->GetTopology()->AtomList.GetAtom(i);
            } else {
                p_root->Atoms[i] = NULL;
            }
        }
        break;
    case O_OR:
        for(int i=0; i < p_root->Owner->GetNumberOfTopologyAtoms(); i++) {
            if( (p_left->Atoms[i] != NULL) || (p_right->Atoms[i] != NULL) ) {
                p_root->Atoms[i] = p_root->Owner->GetTopology()->AtomList.GetAtom(i);
            } else {
                p_root->Atoms[i] = NULL;
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
    case O_AGT:
    case O_ALT: {
        bool result = false;
        switch(p_expr->Modificator) {
        case D_ORIGIN:
            result = p_root->SelectAtomByDistanceFromOrigin(p_expr->Operator,p_expr->Distance);
            break;
        case D_CBOX:
            result = p_root->SelectAtomByDistanceFromCentreOfBox(p_expr->Operator,p_expr->Distance);
            break;
        case D_LIST:
            result = p_root->SelectAtomByDistanceFromList(p_left,p_expr->Operator,p_expr->Distance);
            break;
        case D_COM:
            result = p_root->SelectAtomByDistanceFromCOM(p_left,p_expr->Operator,p_expr->Distance);
            break;
        case D_PLANE:
            result = p_root->SelectAtomByDistanceFromPlane(p_left,p_expr->Operator,p_expr->Distance);
            break;
        default:
            ES_ERROR("not implemented");
            return(false);
        }
        if( result == false ) return(false);
    }
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

bool CAmberMaskASelection::Select(struct SSelection* p_sel)
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
            // no matter of selector - this always means all atoms
            for(int i=0; i < Owner->GetNumberOfTopologyAtoms(); i++) {
                Atoms[i] = Owner->GetTopology()->AtomList.GetAtom(i);
            }
        }
        if( p_item->Index > 0 ) {
            if( p_item->Length >= 1 ) {
                switch(p_sel->Type) {
                case T_RSELECTOR:
                    SelectResidueByIndex(p_item->Index,p_item->Length);
                    break;
                case T_ASELECTOR:
                    SelectAtomByIndex(p_item->Index,p_item->Length);
                    break;
                case T_TSELECTOR:
                    ES_ERROR("index or range for type selector is meaningless");
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
                SelectAtomByName(p_item->Name);
                break;
            case T_TSELECTOR:
                SelectAtomByType(p_item->Name);
                break;
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

void CAmberMaskASelection::SelectAtomByIndex(int index,int length)
{
    for(int i=0; i < length; i++) {
        int aindex = index + i - 1; // transform to index counted from zero

        // check range - this is mandatory because mask can be general
        if( (aindex >= 0) && (aindex < Owner->GetTopology()->AtomList.GetNumberOfAtoms()) ) {
            Atoms[aindex] = Owner->GetTopology()->AtomList.GetAtom(aindex);
        }
    }
}

//------------------------------------------------------------------------------

void CAmberMaskASelection::SelectAtomByName(const char* p_name)
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

    for(int i=0; i < Owner->GetTopology()->AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Owner->GetTopology()->AtomList.GetAtom(i);

        if( strncmp(p_atom->GetName(),p_name,search_len) == 0 ) {
            Atoms[i] = p_atom;
        }

    }
}

//------------------------------------------------------------------------------

void CAmberMaskASelection::SelectAtomByType(const char* p_name)
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

    for(int i=0; i < Owner->GetTopology()->AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Owner->GetTopology()->AtomList.GetAtom(i);
        if( strncmp(p_atom->GetType(),p_name,search_len) == 0 ) {
            Atoms[i] = p_atom;
        }
    }
}

//------------------------------------------------------------------------------

void CAmberMaskASelection::SelectResidueByIndex(int index,int length)
{
    for(int i=0; i < length; i++) {
        int rindex = index + i - 1; // transform to index counted from zero

        // check range - this is mandatory because mask can be general
        if( (rindex >= 0) &&
                (rindex < Owner->GetTopology()->ResidueList.GetNumberOfResidues()) ) {
            CAmberResidue* p_res = Owner->GetTopology()->ResidueList.GetResidue(rindex);

            // select whole residue
            for(int j = 0; j < p_res->GetNumberOfAtoms(); j++) {
                int aindex = j + p_res->GetFirstAtomIndex();
                Atoms[aindex] = Owner->GetTopology()->AtomList.GetAtom(aindex);
            }
        }
    }
}

//------------------------------------------------------------------------------

void CAmberMaskASelection::SelectResidueByName(const char* p_name)
{
    // wild card '*' is treated on parser level!
    // see bool Select(struct SSelection* p_sel) method

    // here we need to manage only '=' wild card
    // this wild card can be only at the end

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
            for(int j = 0; j < p_res->GetNumberOfAtoms(); j++) {
                int aindex = j + p_res->GetFirstAtomIndex();
                Atoms[aindex] = Owner->GetTopology()->AtomList.GetAtom(aindex);
            }
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAmberMaskASelection::SelectAtomByDistanceFromOrigin(SOperator dist_oper,double dist)
{
    if( Owner->GetCoordinates() == NULL ) {
        ES_ERROR("distance operator requires coordinates");
        return(false);
    }

    double dist2 = dist*dist;

    for(int i=0; i < Owner->GetTopology()->AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Owner->GetTopology()->AtomList.GetAtom(i);
        CPoint pos = Owner->GetCoordinates()->GetPosition(i);
        double ldist2 = Square(pos);
        switch(dist_oper) {
        case O_ALT:
            if( ldist2 < dist2 ) Atoms[i] = p_atom;
            break;
        case O_AGT:
            if( ldist2 > dist2 ) Atoms[i] = p_atom;
            break;
        case O_RLT:
        case O_RGT:
        default:
            ES_ERROR("incorrect operator");
            return(false);
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskASelection::SelectAtomByDistanceFromCentreOfBox(SOperator dist_oper,double dist)
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

    for(int i=0; i < Owner->GetTopology()->AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Owner->GetTopology()->AtomList.GetAtom(i);
        CPoint pos = Owner->GetCoordinates()->GetPosition(i);
        double ldist2 = Square(pos-cbox);
        switch(dist_oper) {
        case O_ALT:
            if( ldist2 < dist2 ) Atoms[i] = p_atom;
            break;
        case O_AGT:
            if( ldist2 > dist2 ) Atoms[i] = p_atom;
            break;
        case O_RLT:
        case O_RGT:
        default:
            ES_ERROR("incorrect operator");
            return(false);
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskASelection::SelectAtomByDistanceFromList(CAmberMaskASelection* p_left,
        SOperator dist_oper,double dist)
{
    if( Owner->GetCoordinates() == NULL ) {
        ES_ERROR("distance operator requires coordinates");
        return(false);
    }

    double dist2 = dist*dist;

    for(int i=0; i < Owner->GetTopology()->AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom1 = Owner->GetTopology()->AtomList.GetAtom(i);
        CPoint pos1 = Owner->GetCoordinates()->GetPosition(i);

        for(int j=0; j < Owner->GetTopology()->AtomList.GetNumberOfAtoms(); j++) {
            if( p_left->Atoms[j] == NULL ) continue;
            CPoint pos2 = Owner->GetCoordinates()->GetPosition(j);
            double ldist2 = Square(pos2-pos1);
            bool set = false;
            switch(dist_oper) {
            case O_ALT:
                if( ldist2 < dist2 ) {
                    Atoms[i] = p_atom1;
                    set = true;
                }
                break;
            case O_AGT:
                if( ldist2 > dist2 ) {
                    Atoms[i] = p_atom1;
                    set = true;
                }
                break;
            case O_RLT:
            case O_RGT:
            default:
                ES_ERROR("incorrect operator");
                return(false);
            }
            if( set == true ) break;
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskASelection::SelectAtomByDistanceFromCOM(CAmberMaskASelection* p_left,
        SOperator dist_oper,double dist)
{
    if( Owner->GetCoordinates() == NULL ) {
        ES_ERROR("distance operator requires coordinates");
        return(false);
    }

    // calculate centre of mass
    CPoint     com;
    double     tmass = 0.0;

    for(int i=0; i < Owner->GetTopology()->AtomList.GetNumberOfAtoms(); i++) {
        if( p_left->Atoms[i] == NULL ) continue;
        CAmberAtom* p_atom = Owner->GetTopology()->AtomList.GetAtom(i);
        double mass = p_atom->GetMass();
        com += Owner->GetCoordinates()->GetPosition(i)*mass;
        tmass += mass;
    }

    if( tmass <= 0.0 ) return(true);
    com /= tmass;

    double dist2 = dist*dist;

    // select atoms
    for(int i=0; i < Owner->GetTopology()->AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Owner->GetTopology()->AtomList.GetAtom(i);
        CPoint pos = Owner->GetCoordinates()->GetPosition(i);
        double ldist2 = Square(pos-com);
        switch(dist_oper) {
        case O_ALT:
            if( ldist2 < dist2 ) Atoms[i] = p_atom;
            break;
        case O_AGT:
            if( ldist2 > dist2 ) Atoms[i] = p_atom;
            break;
        case O_RLT:
        case O_RGT:
        default:
            ES_ERROR("incorrect operator");
            return(false);
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskASelection::SelectAtomByDistanceFromPlane(CAmberMaskASelection* p_left,
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

bool CAmberMaskASelection::SelectResidueByDistanceFromOrigin(SOperator dist_oper,double dist)
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
        if( set ) {
            for(int j = 0; j < p_res->GetNumberOfAtoms(); j++) {
                int aindex = j + p_res->GetFirstAtomIndex();
                Atoms[aindex] = Owner->GetTopology()->AtomList.GetAtom(aindex);
            }
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskASelection::SelectResidueByDistanceFromCentreOfBox(SOperator dist_oper,double dist)
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
        if( set ) {
            for(int j = 0; j < p_res->GetNumberOfAtoms(); j++) {
                int aindex = j + p_res->GetFirstAtomIndex();
                Atoms[aindex] = Owner->GetTopology()->AtomList.GetAtom(aindex);
            }
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskASelection::SelectResidueByDistanceFromList(CAmberMaskASelection* p_left,
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
            CPoint pos1 = Owner->GetCoordinates()->GetPosition(aindex);

            for(int k=0; k < Owner->GetTopology()->AtomList.GetNumberOfAtoms(); k++) {
                if( p_left->Atoms[k] == NULL ) continue;
                CPoint pos2 = Owner->GetCoordinates()->GetPosition(k);
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

        // select residue
        if( set ) {
            for(int j = 0; j < p_res->GetNumberOfAtoms(); j++) {
                int aindex = j + p_res->GetFirstAtomIndex();
                Atoms[aindex] = Owner->GetTopology()->AtomList.GetAtom(aindex);
            }
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskASelection::SelectResidueByDistanceFromCOM(CAmberMaskASelection* p_left,
        SOperator dist_oper,double dist)
{
    if( Owner->GetCoordinates() == NULL ) {
        ES_ERROR("distance operator requires coordinates");
        return(false);
    }

    // calculate centre of mass
    CPoint     com;
    double     tmass = 0.0;

    for(int i=0; i < Owner->GetTopology()->AtomList.GetNumberOfAtoms(); i++) {
        if( p_left->Atoms[i] == NULL ) continue;
        CAmberAtom* p_atom = Owner->GetTopology()->AtomList.GetAtom(i);
        double mass = p_atom->GetMass();
        com += Owner->GetCoordinates()->GetPosition(i)*mass;
        tmass += mass;
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
        if( set ) {
            for(int j = 0; j < p_res->GetNumberOfAtoms(); j++) {
                int aindex = j + p_res->GetFirstAtomIndex();
                Atoms[aindex] = Owner->GetTopology()->AtomList.GetAtom(aindex);
            }
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CAmberMaskASelection::SelectResidueByDistanceFromPlane(CAmberMaskASelection* p_left,
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

int CAmberMaskASelection::GetNumberOfSelectedAtoms(void)
{
    if( Atoms == NULL ) return(0);

    int num_of_select_atoms = 0;

    for(int i=0; i < Owner->GetNumberOfTopologyAtoms(); i++) {
        if( Atoms[i] != NULL ) num_of_select_atoms++;
    }

    return(num_of_select_atoms);
}

//------------------------------------------------------------------------------

CAmberAtom* CAmberMaskASelection::GetSelectedAtom(int index)
{
    if( Atoms == NULL ) return(NULL);
    return(Atoms[index]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


