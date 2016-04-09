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
**  Options for NWGOA simulation
*/

#undef NO_HIS
#define HDF5
#define DEFLATE
#define PERFECT_RESTART

/* general */

#define CURVGRID
#define MASKING
#define NONLIN_EOS
#define SOLVE3D
#define SALINITY
#ifdef SOLVE3D
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#endif
#undef FLOATS
#undef STATIONS
#define WET_DRY

/* ice */

#ifdef SOLVE3D
# define  ICE_MODEL
# ifdef ICE_MODEL
#  define ANA_ICE
#  define  ICE_THERMO
#  define  ICE_MK
#  define  ICE_MOMENTUM
#  define  ICE_MOM_BULK
#  define  ICE_EVP
#  define  ICE_STRENGTH_QUAD
#  define  ICE_ADVECT
#  define  ICE_SMOLAR
#  define  ICE_UPWIND
#  define  ICE_BULK_FLUXES
#  define  ICE_I_O
# endif
#endif

/* output stuff */

#define NO_WRITE_GRID
#undef OUT_DOUBLE
#ifndef PERFECT_RESTART
# define RST_SINGLE
#endif
#undef AVERAGES
#undef AVERAGES2
#ifdef SOLVE3D
# undef DIAGNOSTICS_TS
#endif
#undef DIAGNOSTICS_UV

/* advection, dissipation, pressure grad, etc. */

#ifdef SOLVE3D
# define DJ_GRADPS
#endif

#define UV_ADV
#define UV_COR

#ifdef SOLVE3D
# undef TS_A4HADVECTION
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# undef TS_MPDATA
#endif

#define UV_VIS2
#undef UV_SMAGORINSKY
#undef VISC_3DCOEF
#define MIX_S_UV
#define VISC_GRID

#ifdef SOLVE3D
# undef TS_DIF2
# undef MIX_GEO_TS
# undef DIFF_GRID
#endif

/* vertical mixing */

#ifdef SOLVE3D
# define WTYPE_GRID

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

# undef GLS_MIXING
# undef MY25_MIXING

# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  define CRAIG_BANNER
#  define CHARNOK
#  undef GERBI_TKE_FLUX
#  undef AKLIMIT
# endif
#endif

/* surface forcing */

#ifdef SOLVE3D
# define CORE_FORCING
# define BULK_FLUXES
# define CCSM_FLUXES
# if defined BULK_FLUXES || defined CCSM_FLUXES
#  define LONGWAVE_OUT
#  undef DIURNAL_SRFLUX
#  define SOLAR_SOURCE
#  define EMINUSP
#  define ALBEDO_CURVE  /* for water */
#  define ICE_ALB_EC92  /* for ice */
#  undef ALBEDO_CSIM   /* for ice */
#  undef ALBEDO_FILE  /* for both */
#  undef LONGWAVE
# endif
# define SCORRECTION
#else
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
#endif

/* surface and side corrections */

#ifdef SOLVE3D
# define NO_SCORRECTION_ICE
#endif

/* point sources (rivers, line sources) */

/* Not using Runoff now */
#ifdef SOLVE3D
# define ONE_TRACER_SOURCE
#endif

/* tides */

#undef LTIDES
#ifdef LTIDES
# if defined AVERAGES && !defined USE_DEBUG
#  define FILTERED
# endif
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# define TIDES_ASTRO
# undef POT_TIDES
#endif

#define UV_DRAG_GRID
#define ANA_DRAG
#define UV_QDRAG
#define LIMIT_BSTRESS

/* Boundary conditions...careful with grid orientation */

#define RADIATION_2D

/* roms quirks */

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif

