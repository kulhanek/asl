#ifndef AmberBoxH
#define AmberBoxH
/** \ingroup AmberTopology*/
/*! \file AmberBox.hpp */
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
#include <Point.hpp>
#include <ASLMainHeader.hpp>

//---------------------------------------------------------------------------

/// box types

enum EAmberBoxType {
    AMBER_BOX_NONE,
    AMBER_BOX_STANDARD,
    AMBER_BOX_OCTAHEDRAL
};

//---------------------------------------------------------------------------

/// box description for topology

class ASL_PACKAGE CAmberBox {
public:
    CAmberBox(void);
    ~CAmberBox(void);

    /// init all fields - previous data are destroyed
    void InitFields(int iIFBOX);

    /// prepare for new data
    void FreeFields(void);

    /// return type of box
    EAmberBoxType GetType(void);

    /// set the type of box
    void SetType(EAmberBoxType type);

    /// return number of solute atoms
    int GetNumberOfSoluteAtoms(void);

    /// return number of solvent atoms
    int GetNumberOfSolventAtoms(void);

    /// return number of molecules
    int GetNumberOfMolecules(void);

    /// set number of molecules
    bool SetNumberOfMolecules(int number);

    /// set first solvent molecule
    bool SetFirstSolventMolecule(int index);

    /// return number of solute molecules
    int GetNumberOfSoluteMolecules(void);

    /// return number of solvent molecules
    int GetNumberOfSolventMolecules(void);

    /// return number of atoms in particular molecule
    int GetNumberOfAtomsInMolecule(int mol_index);

    /// set number of atoms in particular molecule
    void SetNumberOfAtomsInMolecule(int mol_index,int number);

    /// change list describing number of atoms in each molecule
    /*! method modifies original topology
    */
    bool JoinSoluteMolecules(int number_of_first_molecules);

    /// return the dimensions of box
    CPoint GetBoxDimmensions(void);

    /// set the dimensions of box
    void SetBoxDimmensions(const CPoint& dim);

    /// return the dimensions of box
    CPoint GetBoxAngles(void);

    /// set the dimensions of box
    void SetBoxAngles(const CPoint& angs);

    /// return the box angle
    double GetBoxBeta(void);

    /// set box angle
    void SetBoxBeta(double angle);

    /// update box matrices
    void UpdateBoxMatrices(void);

    /// get box volume
    double GetVolume(void);

    /// get radius of largest inscribed sphere
    double GetLargestSphereRadius(void);

    /// get box center
    const CPoint GetBoxCenter(void);

    /// image point
    const CPoint ImagePoint(const CPoint& pos,
                                         int kx=0,int ky=0,int kz=0,
                                         bool origin=false,
                                         bool familiar=false);

    /// image vector
    const CPoint ImageVector(const CPoint& vec);

    /// overload assigment operator
    void operator = (const CAmberBox& src);

// section of private data ----------------------------------------------------
private:
    /// IFBOX  : set to 1 if standard periodic box, 2 when truncated octahedral
    int     IFBOX;
    int     IPTRES;     //  IPTRES : final residue that is considered
                        //           part of the solute, reset in sander and gibbs
    int     NSPM;       //  NSPM   : total number of molecules
    int     NSPSOL;     //  NSPSOL : the first solvent "molecule"
    int*    NSP;        //  NSP    : the total number of atoms in each
                        //           molecule,necessary to correctly perform the pressure scaling.
    double  DIMM[3];    //  BOX    : the periodic box lengths in the X, Y, and Z directions
    double  ANGS[3];    //  BETA   : periodic box, angle between the XY and YZ planes in degrees.

    double  UCELL[3][3];
    double  RECIP[3][3];

    double  Volume;     // volume of the cell
    double  Radius;     // larges inscribed sphere radius

    int     NumOfSoluteAtoms;

    bool LoadAmber6Box(FILE* p_file);
    bool LoadSolventPointers(FILE* p_file,const char* p_format);
    bool LoadNumsOfMolecules(FILE* p_file,const char* p_format);
    bool LoadBoxInfo(FILE* p_file,const char* p_format);

    bool SaveAmber6Box(FILE* p_file);
    bool SaveSolventPointers(FILE* p_file,const char* p_format);
    bool SaveNumsOfMolecules(FILE* p_file,const char* p_format);
    bool SaveBoxInfo(FILE* p_file,const char* p_format);

    friend class CAmberTopology;
};

//---------------------------------------------------------------------------

#endif
