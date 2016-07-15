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
**  Options for Northeast Pacific (NEP4) simulation
*/
 
/* general */

#define CURVGRID
#define MASKING
#define NONLIN_EOS
#define SOLVE3D
#define SALINITY
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
# undef FLOATS
# undef STATIONS
 
/* ice */

#ifdef SOLVE3D
# define  ICE_MODEL
# define  ICE_THERMO
# define  ICE_MK
# undef   ICE_ALB_EC92
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
#  undef ALBEDO_CLOUD
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
 
/* tides */
 
#undef SSH_TIDES
#undef UV_TIDES
#undef ADD_FSOBC
#undef ADD_M2OBC
#undef RAMP_TIDES
#undef TIDES_ASTRO
#undef POT_TIDES
 
#define RADIATION_2D
 
/* roms quirks */
 
#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif
