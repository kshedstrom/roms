/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group
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
#undef FLOATS
#undef SOLVE3D

#undef UV_ADV
#define UV_COR
#define UV_VIS2
#undef MIX_S_UV
#define  UV_LDRAG

#ifdef SOLVE3D
# undef TS_A4HADVECTION
# undef TS_A4VADVECTION
# undef  TS_U3HADVECTION
# define DJ_GRADPS
# undef TS_DIF2
# undef MIX_S_TS
# define ANA_STFLUX
# define ANA_BTFLUX
# define ANA_VMIX
#endif

#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_FSOBC
#define ANA_M2OBC

#define ANA_PSOURCE
