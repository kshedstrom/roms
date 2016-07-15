/*
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group
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
#undef CCSM_COUPLED
#define CORE_FORCING
#define CCSM_FLUXES
 
#ifdef SOLVE3D
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#endif
#undef FLOATS
#undef STATIONS
#undef WET_DRY

/* ice */

#ifdef SOLVE3D
# undef  ICE_MODEL
#ifdef     ICE_MODEL
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
#endif

/* output stuff */
 
#define NO_WRITE_GRID
#undef OUT_DOUBLE
#define RST_SINGLE
#define AVERAGES
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
#undef UV_SADVECTION
 
#ifdef SOLVE3D
# define TS_U3HADVECTION
# define TS_SVADVECTION
#endif
 
#define UV_VIS2
#define MIX_S_UV
#define VISC_GRID
#undef SPONGE

#ifdef SOLVE3D
# define TS_DIF2
# define MIX_GEO_TS
# define DIFF_GRID
#endif
 
/*#define UV_QDRAG */
#define LIMIT_BSTRESS
 
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
#  define ALBEDO_CURVE
#  undef ALBEDO
#  undef LONGWAVE
# endif
#endif
 
/* surface and side corrections */
 
#ifdef SOLVE3D
# undef SCORRECTION
# undef SRELAXATION /* to do later : correct in set_vbc.F  */
# undef QCORRECTION
#endif

#ifdef SOLVE3D
# undef ANA_TCLIMA
# define ANA_NUDGCOEF
#endif
 
/* point sources (rivers, line sources) */
 
#ifdef SOLVE3D
# undef UV_PSOURCE
# undef TS_PSOURCE
#endif
 
/* tides */
 
#undef LTIDES /* RD test */
#ifdef LTIDES
# define FILTERED
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# undef RAMP_TIDES
# define TIDES_ASTRO
# define POT_TIDES
/*# define UV_LDRAG
# define UV_DRAG_GRID
# define ANA_DRAG
# define LIMIT_BSTRESS
# undef UV_QDRAG */
# define UV_QDRAG 
#else 
# define UV_QDRAG
/* drag could be also something to tweak */ 
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

#define NEMURO
#ifdef NEMURO
# define NEMURO_SAN
# ifdef NEMURO_SAN
#  define FISH_FEEDBACK
#  define PREDATOR
#  define FISHING_FLEET
#  undef  ANA_SPAWN_DIST
#  define EGGS_BISECTION
# endif
# undef  BIO_SEDIMENT
# define HOLLING_GRAZING
# undef  IVLEV_EXPLICIT
# undef  IRON_LIMIT
# undef  IRON_RELAX
# undef  IRON_RSIN
# undef  ANA_BIOLOGY
# define ANA_SPFLUX 
# define ANA_BPFLUX
# undef  ANA_SRFLUX
# undef  DIAGNOSTICS_BIO
#endif

