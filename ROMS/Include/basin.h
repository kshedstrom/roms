/*
** svn $Id: basin.h 975 2009-05-05 22:51:13Z kate $
*******************************************************************************
** Copyright (c) 2002-2009 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Big Bad Basin.
**
** Application flag:   BASIN
** Input script:       ocean_basin.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS4
#define MIX_S_UV
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define SOLVE3D
#define SPLINES
#define EASTERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL
#define BODYFORCE
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
