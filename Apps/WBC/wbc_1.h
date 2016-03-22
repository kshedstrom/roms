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
**  Options for Western Boundary Current example:
*/

#define UV_VIS2
#define UV_COR
#undef AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX

#  if defined WBC_2 || defined WBC_3
#define MASKING
#define ANA_MASK
#  endif

# elif defined CRIT_LAT

/*
**  Options for critical latitude test.
*/

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define SPONGE
#define SOLVE3D
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#undef AVERAGES

#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_M3OBC

#undef  RADIATION_2D
#define NORTH_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_M3RADIATION
#define SOUTH_FSCLAMPED
#define SOUTH_M2CLAMPED
#define SOUTH_M3CLAMPED
