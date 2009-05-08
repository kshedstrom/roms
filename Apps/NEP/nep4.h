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
**  Options for Northeast Pacific (NEP4) simulation
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
 
/* ice */

#ifdef SOLVE3D
# define  ICE_MODEL
# define  ICE_THERMO
# define  ICE_MK
# undef   ICE_ALB_EC92
# undef   ICE_SMOOTH
# define  ICE_MOMENTUM
# define  ICE_MOM_BULK
# define  ICE_EVP
# define  ICE_ADVECT
# define  ICE_SMOLAR
# define  ICE_UPWIND
# define  ICE_BULK_FLUXES
# define  ANA_AIOBC
# define  ANA_HIOBC
# define  ANA_HSNOBC
#endif

/* output stuff */
 
#define NO_WRITE_GRID
#undef OUT_DOUBLE
#define RST_SINGLE
#define AVERAGES
#ifdef SOLVE3D
# undef FILTERED
# define AVERAGES_AKT
# define AVERAGES_AKS
# define AVERAGES_AKV
# define AVERAGES_FLUXES
# undef AVERAGES_QUADRATIC
# undef DIAGNOSTICS_TS
#endif
#undef DIAGNOSTICS_UV
 
/* advection, dissipation, pressure grad, etc. */
 
#ifdef SOLVE3D
# define DJ_GRADPS
#endif
 
#define UV_ADV
#define UV_COR
#define UV_SADVECTION
 
#ifdef SOLVE3D
# define TS_U3HADVECTION
# define TS_SVADVECTION
#endif
 
#define UV_VIS2
#define MIX_S_UV
#define VISC_GRID
#define SPONGE

#ifdef SOLVE3D
# define TS_DIF2
# define MIX_GEO_TS
# define DIFF_GRID
#endif
 
#define UV_QDRAG
 
/* vertical mixing */
 
#ifdef SOLVE3D
# define SOLAR_SOURCE
 
# define LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  define LMD_NONLOCAL
#  define LMD_SHAPIRO
#  undef LMD_DDMIX
#  undef LMD_BKPP
# endif
 
# undef GLS_MIXING
# undef MY25_MIXING
#endif
 
/* surface forcing */
 
#ifdef SOLVE3D
# define BULK_FLUXES
# ifdef BULK_FLUXES
#  define LONGWAVE_OUT
#  define DIURNAL_SRFLUX
#  define EMINUSP
#  undef ANA_SRFLUX
#  undef ALBEDO
#  undef LONGWAVE
# endif
#endif
 
/* surface and side corrections */
 
#ifdef SOLVE3D
# undef SRELAXATION
# undef QCORRECTION
#endif
 
#ifdef SOLVE3D
# undef TCLIMATOLOGY
# undef TCLM_NUDGING
#endif
 
/* point sources (rivers, line sources) */
 
#ifdef SOLVE3D
# define UV_PSOURCE
# define TS_PSOURCE
#endif
 
/* tides */
 
#undef SSH_TIDES
#undef UV_TIDES
#undef ADD_FSOBC
#undef ADD_M2OBC
#undef RAMP_TIDES
#undef TIDES_ASTRO
#undef POT_TIDES
 
/* Boundary conditions...careful with grid orientation */
 
#define EASTERN_WALL
#define NORTHERN_WALL
#undef WESTERN_WALL
#undef SOUTHERN_WALL
 
#define RADIATION_2D
 
#ifndef NORTHERN_WALL
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# ifdef SOLVE3D
#  define NORTH_M3RADIATION
#  define NORTH_M3NUDGING
#  define NORTH_TRADIATION
#  define NORTH_TNUDGING
#  define NORTH_MIGRADIENT
# endif
#endif
 
#ifndef WESTERN_WALL
# define WEST_FSCHAPMAN
# define WEST_M2FLATHER
# ifdef SOLVE3D
#  define WEST_M3RADIATION
#  define WEST_M3NUDGING
#  define WEST_TRADIATION
#  define WEST_TNUDGING
#  define WEST_MIGRADIENT
# endif
#endif
 
#ifndef SOUTHERN_WALL
# define SOUTH_FSCHAPMAN
# define SOUTH_M2FLATHER
# ifdef SOLVE3D
#  define SOUTH_M3RADIATION
#  define SOUTH_M3NUDGING
#  define SOUTH_TRADIATION
#  define SOUTH_TNUDGING
#  define SOUTH_MIGRADIENT
# endif
#endif
 
#ifndef EASTERN_WALL
# define EAST_FSCHAPMAN
# define EAST_M2FLATHER
# ifdef SOLVE3D
#  define EAST_M3RADIATION
#  define EAST_M3NUDGING
#  define EAST_TRADIATION
#  define EAST_TNUDGING
#  define EAST_MIGRADIENT
# endif
#endif
 
/* roms quirks */
 
#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif
