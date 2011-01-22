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
**  Options for South of Java simulation
*/
 
/* general */

#define CURVGRID
#define MASKING
#define NONLIN_EOS
#define SOLVE3D
#define SALINITY
#define SPLINES
# undef FLOATS
# undef STATIONS
 
/* output stuff */
 
#define NO_WRITE_GRID
#undef OUT_DOUBLE
#define RST_SINGLE
#define AVERAGES
#undef DIAGNOSTICS_TS
#undef DIAGNOSTICS_UV
 
/* advection, dissipation, pressure grad, etc. */
 
#define DJ_GRADPS
 
#define UV_ADV
#define UV_COR
#define UV_SADVECTION
 
#define TS_U3HADVECTION
#define TS_SVADVECTION
 
#define UV_VIS2
#define MIX_S_UV
#define VISC_GRID

          
#define TS_DIF2
#define MIX_GEO_TS
#define DIFF_GRID
 
#define UV_QDRAG
 
/* vertical mixing */
 
#define SOLAR_SOURCE
 
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#  undef LMD_DDMIX
#  undef LMD_BKPP
#  undef LMD_SHAPIRO
#endif
 
# undef GLS_MIXING
# undef MY25_MIXING
 
/* surface forcing */
 
#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE_OUT
# define DIURNAL_SRFLUX
# define EMINUSP
# undef ANA_SRFLUX
# undef ALBEDO
# undef LONGWAVE
#endif
 
/* surface and side corrections */
 
# undef SRELAXATION
# undef QCORRECTION
 
#define SPONGE
# undef TCLIMATOLOGY
# undef TCLM_NUDGING
 
/* Boundary conditions...careful with grid orientation */
 
#undef EASTERN_WALL
#define NORTHERN_WALL
#undef WESTERN_WALL
#undef SOUTHERN_WALL
 
#define RADIATION_2D
 
#ifndef NORTHERN_WALL
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# define NORTH_M3RADIATION
# define NORTH_M3NUDGING
# define NORTH_TRADIATION
# define NORTH_TNUDGING
#endif
 
#ifndef WESTERN_WALL
# define WEST_FSCHAPMAN
# define WEST_M2FLATHER
# define WEST_M3RADIATION
# define WEST_M3NUDGING
# define WEST_TRADIATION
# define WEST_TNUDGING
#endif
 
#ifndef SOUTHERN_WALL
# define SOUTH_FSCHAPMAN
# define SOUTH_M2FLATHER
# define SOUTH_M3RADIATION
# define SOUTH_M3NUDGING
# define SOUTH_TRADIATION
# define SOUTH_TNUDGING
#endif
 
#ifndef EASTERN_WALL
# define EAST_FSCHAPMAN
# define EAST_M2FLATHER
# define EAST_M3RADIATION
# define EAST_M3NUDGING
# define EAST_TRADIATION
# define EAST_TNUDGING
#endif
 
/* roms quirks */
 
#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif
