#ifndef AmberRestartH
#define AmberRestartH
/** \ingroup AmberRestart*/
/*! \file AmberRestart.hpp */
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

#include <ASLMainHeader.hpp>
#include <stdio.h>
#include <Point.hpp>
#include <SmallString.hpp>

//---------------------------------------------------------------------------

class CAmberTopology;
class CXMLElement;

//---------------------------------------------------------------------------

/// restart file IO class

class ASL_PACKAGE CAmberRestart {
public:
    CAmberRestart(void);
    ~CAmberRestart(void);

    /// assign topology with restart file
    /*! previous data are discarded !
    */
    void AssignTopology(CAmberTopology* p_topology);

    /// create and clear all data
    bool Create(void);

    /// relase all data
    void Release(void);

    /// get number of atoms
    int GetNumberOfAtoms(void) const;

    /// get topology
    CAmberTopology* GetTopology(void);

    /// load coordinates - velocities and box if present
    bool Load(const CSmallString& name,bool allow_stdin=false);

    /// save coordinates - velocities and box if present
    bool Save(const CSmallString& name,bool allow_stdout=false);

    /// load coordinates - velocities and box if present
    bool Load(FILE *fin);

    /// save coordinates - velocities and box if present
    bool Save(FILE *fout);

    /// load only coordinates from XML record
    bool LoadSnapshot(CXMLElement* p_ele);

    /// load only velocities from XML record (CRD->VEL)
    bool LoadVelocities(CXMLElement* p_ele);

    /// save only coordinates to XML record
    void SaveSnapshot(CXMLElement* p_ele);

    /// return position of atom of index
    const CPoint& GetPosition(int index) const;

    /// return mass of atom of index
    double GetMass(int index) const;

    /// set position of atom of index
    void SetPosition(int index,const CPoint& pos);

    /// return velocity of atom of index
    const CPoint& GetVelocity(int index) const;

    /// set velocity of atom of index
    void SetVelocity(int index,const CPoint& vel);

    /// return box size
    const CPoint& GetBox(void) const;

    /// set new box size
    void SetBox(const CPoint& newbox);

    /// return box angles
    const CPoint& GetAngles(void) const;

    /// set box angles
    void SetAngles(const CPoint& newangles);

    /// return true if velocities were loaded from restart file
    bool AreVelocitiesLoaded(void) const;

    /// return true if box is present
    bool IsBoxPresent(void) const;

    /// get snapshot time
    double GetTime(void);

    /// get snapshot time
    void SetTime(const double time);

    /// get title
    const CSmallString& GetTitle(void);

    /// set title
    void SetTitle(const CSmallString& title);

    /// overload assigment operator
    void operator = (const CAmberRestart& src);

// section of private data ----------------------------------------------------
private:
    CAmberTopology* Topology;
    CSmallString    Title;
    double          Time;
    int             NumberOfAtoms;
    CPoint*         Positions;
    CPoint*         Velocities;
    CPoint          Box;    // box dimensions
    CPoint          Box1;   // box angles
    bool            VelocitiesLoaded;

    static CPoint zero;

    friend class CAmberTrajectory;
};
//---------------------------------------------------------------------------

#endif
