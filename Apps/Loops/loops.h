/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2011 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for overturning loopy flow
*/

#undef STATIONS
#undef UV_ADV
#undef UV_COR
#undef TS_A4HADVECTION
#undef TS_A4VADVECTION
#undef  TS_U3HADVECTION
#define DJ_GRADPS
#undef UV_VIS2
#undef TS_DIF2
#undef MIX_S_UV
#undef MIX_S_TS
#define SOLVE3D
#define FLOATS
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_VMIX
#define EASTERN_WALL
#define WESTERN_WALL
#define NS_PERIODIC
#define  UV_LDRAG

