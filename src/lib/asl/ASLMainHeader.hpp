#ifndef ASLMainHeaderH
#define ASLMainHeaderH
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

#include <HiPoLyMainHeader.hpp>

//------------------------------------------------------------------------------

#define ASL_VERSION "ASL 2.0.SVNVERSION (DATE)"

//------------------------------------------------------------------------------

extern const char* LibBuildVersion_ASL;

//------------------------------------------------------------------------------

#if defined _WIN32 || defined __CYGWIN__
#ifdef ASL_BUILDING_DLL
#ifdef __GNUC__
#define ASL_DLL_PUBLIC __attribute__((dllexport))
#else
#define ASL_DLL_PUBLIC __declspec(dllexport)
#endif
#else
#ifdef __GNUC__
#define ASL_DLL_PUBLIC __attribute__((dllimport))
#else
#define ASL_DLL_PUBLIC __declspec(dllimport)
#endif
#define ASL_DLL_LOCAL
#endif
#else
#if __GNUC__ >= 4
#define ASL_DLL_PUBLIC __attribute__ ((visibility("default")))
#define ASL_DLL_LOCAL  __attribute__ ((visibility("hidden")))
#else
#define ASL_DLL_PUBLIC
#define ASL_DLL_LOCAL
#endif
#endif

#define ASL_PACKAGE ASL_DLL_PUBLIC

//---------------------------------------------------------------------------

#endif
