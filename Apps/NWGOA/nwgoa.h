/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2015 The ROMS/TOMS Group
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
#undef HDF5
#undef DEFLATE
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
/* FOR TIDES TEST # define  ICE_MODEL */
# ifdef ICE_MODEL
#  undef ANA_ICE
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
#undef SPONGE

#ifdef SOLVE3D
# undef TS_DIF2
# undef MIX_GEO_TS
# undef DIFF_GRID
#endif

/* vertical mixing */

#ifdef SOLVE3D
# define WTYPE_GRID

# undef LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  undef LMD_BKPP
#  define LMD_NONLOCAL
#  define LMD_SHAPIRO
#  undef LMD_DDMIX
# endif

# define GLS_MIXING
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

/* #ifdef SOLVE3D */
# ifdef FOOOO /* FOR TIDES TEST */
# define CORE_FORCING
# define BULK_FLUXES
# define CCSM_FLUXES
# if defined BULK_FLUXES || defined CCSM_FLUXES
#  define LONGWAVE_OUT
#  undef DIURNAL_SRFLUX
#  define SOLAR_SOURCE
#  define EMINUSP
#  undef ANA_SRFLUX
#  undef ALBEDO_CLOUD
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

/* Using Runoff now */
#ifdef SOLVE3D
# undef RUNOFF
# define ONE_TRACER_SOURCE
#endif

/* tides */

#define LTIDES
#ifdef LTIDES
# define FILTERED
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# undef RAMP_TIDES
# define TIDES_ASTRO
# undef POT_TIDES

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
#undef BIO_UMAINE

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

#ifdef BIO_UMAINE
# define CARBON
# define OXYGEN
# define PRIMARY_PROD
# define SINK_OP2
# define TALK_NONCONSERV
# undef OPTIC_UMaine
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
# define IRON_LIMIT        /* Add iron as passive Nth tracer */
# define IRON_RELAX
# undef  IRON_RSIN
#endif
