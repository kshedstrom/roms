/*
** svn $Id: windbasin.h 975 2009-05-05 22:51:13Z kate $
*******************************************************************************
** Copyright (c) 2002-2009 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Wind-Driven Constant Coriolis Basin Test.
**
** Application flag:   WINDBASIN
** Input script:       ocean_windbasin.in
*/

#undef UV_ADV
#define UV_COR
#define UV_QDRAG
#define SOLVE3D
#define SPLINES
#define AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX

