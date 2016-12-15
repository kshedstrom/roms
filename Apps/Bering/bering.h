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

#undef ROMS_MODEL
#undef WRF_MODEL
#undef MCT_LIB
#undef ATM2OCN_FLUXES  /* not sure about this with ice */
#undef NO_LBC_ATT

#define NO_HIS
#undef HDF5
#undef DEFLATE
#undef PARALLEL_IN
#undef PARALLEL_OUT
#define PERFECT_RESTART

/* general */

#define CURVGRID
#define MASKING
#define SOLVE3D
#ifdef SOLVE3D
# define SALINITY
# define NONLIN_EOS
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#endif
#undef FLOATS
#define STATIONS
#define WET_DRY

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
# define CICE_MODEL
# ifdef CICE_MODEL
#  define SNOWFALL
#  define SNOW_FROM_RAIN
# endif

# undef  ICE_MODEL
# ifdef ICE_MODEL
#  define ANA_ICE
#  undef  OUTFLOW_MASK
#  undef  FASTICE_CLIMATOLOGY
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
#  define ICE_CONVSNOW
#  define ICE_I_O
#  undef  MELT_PONDS
#  undef  ANA_AIOBC
#  undef  ANA_HIOBC
#  undef  ANA_HSNOBC
# endif
#endif

/* output stuff */

#define NO_WRITE_GRID
#undef OUT_DOUBLE
#ifndef PERFECT_RESTART
# define RST_SINGLE
#endif
#define AVERAGES
#undef AVERAGES2
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
#undef UV_SADVECTION

#ifdef SOLVE3D
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# undef TS_MPDATA
#endif

#define UV_VIS2
#undef VISC_3DCOEF
#define MIX_S_UV
#define VISC_GRID
#undef SPONGE

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

# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  define AKLIMIT
# endif
#endif

/* surface forcing */

#ifdef SOLVE3D
# ifndef ATM2OCN_FLUXES
#  define CORE_FORCING
#  define BULK_FLUXES
#  define CCSM_FLUXES
#  undef ARCTIC_MERRA_HACK
# endif
# if defined BULK_FLUXES || defined CCSM_FLUXES
#  define LONGWAVE_OUT
#  undef DIURNAL_SRFLUX
#  define SOLAR_SOURCE
#  define EMINUSP
#  undef ANA_SRFLUX
#  undef ALBEDO_CLOUD
#  define ALBEDO_CURVE  /* for water */
#  undef ICE_ALB_EC92  /* for ice */
#  define ALBEDO_CSIM   /* for ice */
#  undef ALBEDO_FILE  /* for both */
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

/* Using Runoff now */
#ifdef SOLVE3D
# define RUNOFF
#endif

/* tides */

#define LTIDES
#ifdef LTIDES
# if defined AVERAGES && !defined USE_DEBUG
#  define FILTERED
# endif
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# undef RAMP_TIDES
# define TIDES_ASTRO
# define POT_TIDES

# undef UV_LDRAG
# define UV_DRAG_GRID
# define ANA_DRAG
# define LIMIT_BSTRESS
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

/*
**  Biological model options.
*/
#undef NEMURO

#if defined NEMURO
# define BIO_SEDIMENT
# define NEMURO_SED1
# undef ANA_BIOLOGY       /* analytical biology initial conditions */
# define IRON_LIMIT        /* Add iron as passive 11th tracer */
# define IRON_RELAX
# undef  IRON_RSIN
# define HOLLING_GRAZING
# undef  IVLEV_EXPLICIT
# undef  ANA_BIOSWRAD
# undef  DIAGNOSTICS_BIO
#endif
