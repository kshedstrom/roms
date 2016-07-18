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
**  Options for ARCTIC simulation
*/

#undef NO_HIS
#define GLOBAL_PERIODIC
#define HDF5
#define DEFLATE
#undef PARALLEL_IN
#undef PARALLEL_OUT
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
#define FLOATS
#define STATIONS
#undef WET_DRY

#undef T_PASSIVE
#ifdef T_PASSIVE
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX 
# define ANA_PASSIVE
# define TRC_PSOURCE
# define ANA_TRC_PSOURCE
# define AGE_MEAN
#endif

/* ice */

#ifdef SOLVE3D
# undef CICE_MODEL
# ifdef CICE_MODEL
#  define SNOWFALL
#  define SNOW_FROM_RAIN
#  define INI_GLORYS_ICE
#  undef ICE_LOG_LAYER
# endif

# define  ICE_MODEL
# ifdef ICE_MODEL
#  define ANA_ICE
#  define INI_GLORYS_ICE
#  define  OUTFLOW_MASK
#  undef  FASTICE_CLIMATOLOGY
#  define ICE_LANDFAST
#  define ICE_THERMO
#  define ICE_MK
#  define ICE_MOMENTUM
#  define ICE_MOM_BULK
#  define ICE_EVP
#  define ICE_STRENGTH_QUAD
#  define ICE_ADVECT
#  define ICE_SMOLAR
#  define ICE_UPWIND
#  define ICE_BULK_FLUXES
#  define ICE_CONVSNOW
#  define ICE_I_O
# endif
#endif

/* output stuff */

#define NO_WRITE_GRID
#undef OUT_DOUBLE
#ifndef PERFECT_RESTART
# define RST_SINGLE
#endif
#define AVERAGES
#define AVERAGES2
#ifdef SOLVE3D
# undef AVERAGES_DETIDE
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
# define TS_DIF2
# define MIX_GEO_TS
# define DIFF_GRID
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
# endif
#endif

/* surface forcing */

#ifdef SOLVE3D
# define CORE_FORCING
# define BULK_FLUXES
# define CCSM_FLUXES
# undef ARCTIC_MERRA_HACK
# if defined BULK_FLUXES || defined CCSM_FLUXES
#  define LONGWAVE_OUT
#  undef DIURNAL_SRFLUX
#  define SOLAR_SOURCE
#  define EMINUSP
#  undef ALBEDO_CLOUD
#  define ALBEDO_CURVE  /* for water */
#  undef ICE_ALB_EC92   /* for ice */
#  define ALBEDO_CSIM   /* for ice */
#  undef ALBEDO_FILE    /* for both */
#  undef LONGWAVE
# endif
#endif

/* surface and side corrections */

#ifdef SOLVE3D
# define SCORRECTION
# define NO_SCORRECTION_ICE
# undef QCORRECTION
#endif

#ifdef SOLVE3D
# define ANA_NUDGCOEF
#endif

/* point sources (rivers, line sources) */

/* Using both for different regions */
#ifdef SOLVE3D
# define RUNOFF
# define TWO_D_TRACER_SOURCE
#endif

/* tides */

#define LTIDES
#ifdef LTIDES
# if defined AVERAGES && !defined USE_DEBUG
#  undef FILTERED
# endif
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# undef RAMP_TIDES
# define TIDES_ASTRO
# define POT_TIDES

# define UV_DRAG_GRID
# define ANA_DRAG
# define LIMIT_BSTRESS
# undef UV_LDRAG
# define UV_QDRAG
#else
# define UV_QDRAG
#endif

/* Boundary conditions...careful with grid orientation */

#define RADIATION_2D

/* roms quirks */

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif
