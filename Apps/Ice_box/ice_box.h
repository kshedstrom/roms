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
**  Options for BERING simulation
*/

/* COAWST stuff */
#undef ROMS_MODEL
#undef WRF_MODEL
#undef MCT_LIB
#undef ATM2OCN_FLUXES  /* not sure about this with ice */
#undef NO_LBC_ATT

/* ROMS stuff */
#define NO_HIS
#undef HDF5
#undef DEFLATE
#define PERFECT_RESTART

/* general */

#define SOLVE3D
#ifdef SOLVE3D
# define SALINITY
# define NONLIN_EOS
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
# undef TS_FIXED
#endif
#define STATIONS
#define ANA_INITIAL
#define ANA_GRID
#define SPHERICAL

/* mixing */

# define LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  undef LMD_BKPP
#  define LMD_NONLOCAL
#  define LMD_SHAPIRO
#  undef LMD_DDMIX
# endif

/* ice */

#ifdef SOLVE3D
# undef CICE_MODEL
# ifdef CICE_MODEL
#  define SNOWFALL
#  define SNOW_FROM_RAIN
# endif

# define  ICE_MODEL
# undef NO_SNOW
# ifdef ICE_MODEL
#  define ANA_ICE
#  define SNOWFALL
#  define  ICE_THERMO
#  define  ICE_MK
#  undef  ICE_MOMENTUM
#  undef  ICE_MOM_BULK
#  undef  ICE_EVP
#  undef  ICE_STRENGTH_QUAD
#  define  ICE_ADVECT   /* Note that we need these two for the */
#  define  ICE_SMOLAR   /* timestepping to work correctly.     */
#  undef  ICE_UPWIND
#  define  ICE_BULK_FLUXES
#  define ICE_CONVSNOW
#  undef  MELT_PONDS
#  define ICE_I_O
# endif
#endif

/* output stuff */

#define NO_WRITE_GRID
#undef OUT_DOUBLE
#ifndef PERFECT_RESTART
# define RST_SINGLE
#endif
#undef AVERAGES

/* advection, dissipation, pressure grad, etc. */

#ifdef SOLVE3D
# define DJ_GRADPS
#endif

#undef UV_ADV
#undef UV_COR
#undef UV_SADVECTION

#ifdef SOLVE3D
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# undef TS_MPDATA
#endif

#undef UV_VIS2
#undef MIX_S_UV

/* surface forcing */

#ifdef SOLVE3D
# ifndef ATM2OCN_FLUXES
#  define CORE_FORCING
#  define BULK_FLUXES
#  define CCSM_FLUXES
#  undef ARCTIC_MERRA_HACK
# endif
# if defined BULK_FLUXES
#  define LONGWAVE_OUT
#  define SOLAR_SOURCE
#  define EMINUSP
#  define ANA_LRFLUX
#  define ANA_SRFLUX
#  undef ANA_ALBEDO
#  define ALBEDO
#  define ANA_SNOW
#  undef ALBEDO_CLOUD
#  undef ALBEDO_CURVE  /* for water */
#  undef ICE_ALB_EC92  /* for ice */
#  undef ALBEDO_CSIM   /* for ice */
#  undef ALBEDO_FILE  /* for both */
#  undef LONGWAVE
#  define ANA_PAIR
#  define ANA_TAIR
#  define ANA_HUMIDITY
#  define ANA_WINDS
#  define ANA_CLOUD
#  define ANA_RAIN
# endif
#endif

#define UV_QDRAG

/* roms quirks */

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif

