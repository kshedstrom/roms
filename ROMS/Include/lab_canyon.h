/*
** svn $Id: lab_canyon.h 975 2009-05-05 22:51:13Z kate $
*******************************************************************************
** Copyright (c) 2002-2009 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Lab Canyon (Polar Coordinates).
**
** Application flag:   LAB_CANYON
** Input script:       ocean_lab_canyon.in
*/

#define UV_COR
#define UV_ADV
#define UV_LDRAG
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define SOLVE3D
#define CURVGRID
#define AVERAGES
#define SPLINES
#define NS_PERIODIC
#define EASTERN_WALL
#define WESTERN_WALL
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
