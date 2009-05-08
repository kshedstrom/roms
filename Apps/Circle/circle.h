/*
** svn $Id: basin.h 8 2007-02-06 19:00:29Z arango $
*******************************************************************************
** Copyright (c) 2002-2009 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for Circle:
*/

#define UV_VIS2
#define UV_COR
#undef UV_QDRAG
#define UV_LDRAG

#undef AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define MASKING
#define ANA_MASK
