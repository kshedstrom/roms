/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2013 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for 3 km Coastal Gulf of Alaska
*/

#define SOLVE3D
#ifdef SOLVE3D
# if defined BERING
#  define  ICE_MODEL
#  undef  ICE_THERMO
#  undef  ICE_MK
#  undef   ICE_ALB_EC92
#  undef   ICE_SMOOTH
#  define  ICE_MOMENTUM
#  define  ICE_MOM_BULK
#  define  ICE_EVP
#  define  ICE_ADVECT
#  define  ICE_SMOLAR
#  define  ICE_UPWIND
!#  define  ICE_BULK_FLUXES
!#  define  BULK_FLUXES
#  define  CLOUDS
# else
#  define SCORRECTION
#  define QCORRECTION
#  define BULK_FLUXES
#  define LONGWAVE
#  define ALBEDO
#  define ANA_SRFLUX
# endif
#else
# define ANA_INITIAL
#endif

#undef OFFLINE_FLOATS

#define SPLINES
#undef FLOATS
#define MASKING
#define UV_ADV
#define UV_COR
#define UV_QDRAG
#ifdef SOLVE3D
# define DJ_GRADPS
# define TS_U3HADVECTION
# define TS_SVADVECTION
# undef  TS_C2HADVECTION
# undef  TS_C2VADVECTION
#endif

!#if defined OFFLINE_FLOATS
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
!# define ANA_VMIX
# define ANA_SSFLUX
!#else
# ifdef SOLVE3D
#  define ADD_FSOBC
#  define ADD_M2OBC
#  define LMD_MIXING
#  define LMD_RIMIX
#  undef LMD_SHAPIRO
#  undef LMD_CONVEC
#  define LMD_SKPP
#  define LMD_BKPP
#  define LMD_NONLOCAL
#  undef MY25_MIXING
#  define SOLAR_SOURCE
#  undef NCEP_FORCE
#  if defined NCEP_FORCE || defined BULK_FLUXES
#   define LONGWAVE
#   define ALBEDO
#   define ANA_SRFLUX
#  else
#   define SCORRECTION
#   define QCORRECTION
#  endif
# endif
# define SSH_TIDES
# define UV_TIDES
# undef TIDES_ASTRO
!#endif

#define UV_VIS2
#define MIX_S_UV
#define  STATIONS        /* define if writing out station data */
#undef  STATIONS_CGRID  /* define if extracting data at native C-grid */

#define NO_WRITE_GRID
#define NONLIN_EOS
#define CURVGRID
#undef AVERAGES
#ifdef SOLVE3D
# define SPONGE
#endif

#ifdef SOLVE3D
# define SALINITY
# undef TS_DIF2
# undef MIX_GEO_TS
# define FILTERED

# undef TCLIMATOLOGY
# undef TCLM_NUDGING

#endif

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX     /* analytical surface momentum stress */
# define ANA_INITIAL
#endif

#ifndef OFFLINE_FLOATS
# define RADIATION_2D
#endif
