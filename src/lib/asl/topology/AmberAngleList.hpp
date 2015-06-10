#ifndef AmberAngleListH
#define AmberAngleListH
/** \ingroup AmberTopology*/
/*! \file AmberAngleList.hpp */
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
#include <ASLMainHeader.hpp>
#include <AmberAngleType.hpp>
#include <AmberAngle.hpp>

//---------------------------------------------------------------------------

/// list of angles for topology

class ASL_PACKAGE CAmberAngleList {
public:
    CAmberAngleList(void);
    ~CAmberAngleList(void);

    /// init all fields - old data are destroyed
    void InitFields(int iNTHETH, int iMTHETA, int iNUMANG,
                                 int iNTHETA, int iNGPER, int iMGPER);

    /// prepare for new data
    void FreeFields(void);

    /// return number of angles not perturbed
    /*! AMBER abreviation: NTHETH + MTHETA
    */
    int GetNumberOfAngles(void) const;

    /// return number of angles containing hydrogen and not perturbed
    /*! AMBER abreviation: NTHETH
    */
    int GetNumberOfAnglesWithHydrogen(void) const;

    /// return number of angles not containing hydrogen and not perturbed
    /*! AMBER abreviation: MTHETA
    */
    int GetNumberOfAnglesWithoutHydrogen(void) const;

    /// return number of perturbed angles
    /*! AMBER abreviation: NGPER
    */
    int GetNumberOfPerturbedAngles(void) const;

    /// return number of completely perturbed angles
    /*! AMBER abreviation: MGPER
    */
    int GetNumberOfCompletelyPerturbedAngles(void) const;

    /// return number of angle types
    /*! AMBER abreviation: NUMANG
    */
    int GetNumberOfAngleTypes(void) const;

// no longer supported parameters ------------------------------------------
    /// return NTHETA value
    int GetNTHETA(void) const;

// pointers ----------------------------------------------------------------
    /// return pointer to angle containing hydrogen
    CAmberAngle*     GetAngleWithHydrogen(int index) const;

    /// return pointer to angle not containing hydrogen
    CAmberAngle*     GetAngleWithoutHydrogen(int index) const;

    /// return pointer to angle
    /*!
        method return pointer to any non-perturbed angle
        index from
            0 to GetNumberOfAnglesWithoutHydrogen() - 1 points
            to angle not containing hydrogen
        index from
            GetNumberOfAnglesWithoutHydrogen() to GetNumberOfAngles() - 1 points
            to angle containing hydrogen
    */
    CAmberAngle*     GetAngle(int index) const;

    /// return pointer to perturbed angle
    CAmberAngle*     GetPerturbedAngle(int index) const;

    /// return pointer to angle type
    CAmberAngleType* GetAngleType(int index) const;

    /// overload assigment operator
    void operator = (const CAmberAngleList& src);

// private methods ----------------------------------------------------------
private:
    int NTHETH;   //NTHETH : number of angles containing hydrogen
    int MTHETA;   //MTHETA : number of angles not containing hydrogen
    int NUMANG;   //NUMANG : number of unique angle types
    int NTHETA;   //NTHETA : MTHETA + number of constraint angles
    int NGPER;    //NGPER  : number of angles to be perturbed
    int MGPER;    //MGPER  : number of angles with atoms completely in perturbed group

    CAmberAngleType* AngleTypes;
    CAmberAngle*     AngleWithoutHydrogens;
    CAmberAngle*     AngleWithHydrogens;
    CAmberAngle*     PerturbedAngles;

    bool LoadAngleTK(FILE* p_file,const char* p_format);
    bool LoadAngleTEQ(FILE* p_file,const char* p_format);
    bool LoadAnglesWithHydrogens(FILE* p_file,const char* p_format);
    bool LoadAnglesWithoutHydrogens(FILE* p_file,const char* p_format);
    bool LoadPerturbedAngles(FILE* p_file,const char* p_format);
    bool LoadPerturbedAngleTypeIndexes(FILE* p_file,const char* p_format);

    bool SaveAngleTK(FILE* p_file,const char* p_format);
    bool SaveAngleTEQ(FILE* p_file,const char* p_format);
    bool SaveAnglesWithHydrogens(FILE* p_file,const char* p_format);
    bool SaveAnglesWithoutHydrogens(FILE* p_file,const char* p_format);
    bool SavePerturbedAngles(FILE* p_file,const char* p_format);
    bool SavePerturbedAngleTypeIndexes(FILE* p_file,const char* p_format);

    friend class CAmberTopology;
};

//---------------------------------------------------------------------------

#endif
