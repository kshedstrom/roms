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
**  Options for Northeast Pacific (NEP5) simulation
*/

#define NO_HIS
#undef NETCDF4
#undef PARALLEL_IO
#undef OFFLINE_FLOATS

/* general */

#define CURVGRID
#define MASKING
#define NONLIN_EOS
#define SOLVE3D
#define SALINITY
#ifdef SOLVE3D
# undef SPLINES_VDIFF
# undef SPLINES_VVISC
# undef RI_SPLINES
#endif
#define FLOATS
#define STATIONS
#undef WET_DRY

#undef T_PASSIVE
#ifdef T_PASSIVE
# define ANA_PASSIVE
# define TRC_PSOURCE
# define ANA_TRC_PSOURCE
# define AGE_MEAN
#endif

/* ice */

#ifdef SOLVE3D
# define  ICE_MODEL
# ifdef ICE_MODEL
#  define  ICE_THERMO
#  define  ICE_MK
#  undef   ICE_ALB_EC92
#  define  ICE_MOMENTUM
#  define  ICE_MOM_BULK
#  define  ICE_EVP
#  define  ICE_ADVECT
#  define  ICE_SMOLAR
#  define  ICE_UPWIND
#  define  ICE_BULK_FLUXES
#  define  ANA_AIOBC
#  define  ANA_HIOBC
#  define  ANA_HSNOBC
# endif
#endif

/* output stuff */

#define NO_WRITE_GRID
#undef OUT_DOUBLE
#define RST_SINGLE
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
#define UV_SADVECTION

#ifdef SOLVE3D
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# undef TS_MPDATA
#endif

#define UV_VIS2
#undef UV_SMAGORINSKY
#define VISC_3DCOEF
#define MIX_S_UV
#define VISC_GRID
#undef SPONGE

#ifdef SOLVE3D
# define TS_DIF2
# define MIX_GEO_TS
# define DIFF_GRID
#endif

/* vertical mixing */

#ifdef SOLVE3D
# define SOLAR_SOURCE
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
# if defined BULK_FLUXES || defined CCSM_FLUXES
#  define LONGWAVE_OUT
#  define DIURNAL_SRFLUX
#  define EMINUSP
#  undef ANA_SRFLUX
#  undef ALBEDO_CLOUD
#  define ALBEDO_CURVE
#  undef LONGWAVE
# endif
#endif

/* surface and side corrections */

#ifdef SOLVE3D
# define SCORRECTION
# undef SRELAXATION
# undef QCORRECTION
#endif

#ifdef SOLVE3D
# undef TCLIMATOLOGY
# undef TCLM_NUDGING
#endif

/* point sources (rivers, line sources) */

/* Using Runoff instead now */
#ifdef SOLVE3D
# define RUNOFF
# undef ANA_PSOURCE
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
# define POT_TIDES

# define UV_LDRAG
# define UV_DRAG_GRID
# define ANA_DRAG
# define LIMIT_BSTRESS
# undef UV_QDRAG
#else
# define UV_QDRAG
#endif

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
#undef BIO_GOANPZ        /* Sarah Hinckley's 11 box model */
#undef BEST_NPZ         /* Georgina Gibsons BEST NPZ model  */

#if defined BEST_NPZ || defined BIO_GOANPZ
# undef  BIOFLUX           /* sum Nitrogen fluxes between boxes */
# define ANA_BIOLOGY       /* analytical biology initial conditions */
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
# define DIAPAUSE          /* Enable Neocalanus seasonal vertical migration */
# undef FLOAT_VWALK
# define IRON_LIMIT        /* Add iron as passive 11th tracer */
# undef TCLM_NUDGING      /* Nudging of tracer climatology for iron */
#endif

#if defined NEMURO
# define BIO_SEDIMENT
# define NEMURO_SED1
# undef ANA_BIOLOGY       /* analytical biology initial conditions */
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
# define IRON_LIMIT        /* Add iron as passive 11th tracer */
# define IRON_RELAX
# undef  IRON_RSIN
# define HOLLING_GRAZING
# undef  IVLEV_EXPLICIT
# undef  ANA_BIOSWRAD
# undef  DIAGNOSTICS_BIO
#endif

#ifdef BEST_NPZ
# define        NEWSHADE    /* Use Craig''s formulation for self shading in PAR calc+                       Else use Sarah''s self-shading from original NPZ code */
# undef        KODIAK_IRAD /* Generate irradiance with curve matching Kodiak data
                       Else use shortwave radiation (srflx) as irradiance   */
# define JELLY
# define STATIONARY
# define STATIONARY2
# define PROD3
# define PROD2
# define BENTHIC /*FENNEL or BENTHIC or TRAP*/
# define ICE_BIO
# undef CLIM_ICE_1D

# undef SINKVAR      /* for variable sinking rate*/
# undef DENMAN

# undef OFFLINE_BIOLOGY   /* define if offline simulation of bio tracers */
#   if defined OFFLINE_BIOLOGY
#    define AKSCLIMATOLOGY   /* Processing of AKS climatology */
#    undef ANA_AKSCLIMA      /* Processing of AKS climatology */
#   endif
#  undef DIAPAUSE          /* Enable Neocalanus seasonal vertical migration */
#  define  IRON_LIMIT        /* Add iron as passive 13th tracer */
#    if defined IRON_LIMIT || defined CLIM_ICE_1D
#      if !defined OFFLINE_BIOLOGY
#       define TCLM_NUDGING    /* Nudging of tracer climatology for iron */
#       undef  ANA_TCLIMA     /* analytical tracers climatology for iron */
#       define TCLIMATOLOGY   /* Processing of tracer climatology for iron */
#      endif
#    endif
#endif

