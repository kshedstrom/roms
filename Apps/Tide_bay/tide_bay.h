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
**  Tidal half pipe.
*/

#undef STATIONS
#define FLOATS
#undef  DIAGNOSTICS_UV

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define UV_LDRAG

#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX

#define ANA_FSOBC
#define ANA_M2OBC
