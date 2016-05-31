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
#undef SOLVE3D
#ifdef SOLVE3D
# define NONLIN_EOS
# define SALINITY
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#endif
#undef FLOATS
#define STATIONS
#define WET_DRY

#undef T_PASSIVE
#ifdef T_PASSIVE
# define ANA_PASSIVE
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
#undef RST_SINGLE
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
#define UV_VIS2

#ifdef SOLVE3D
# undef UV_SMAGORINSKY
# define VISC_3DCOEF
# define MIX_S_UV
# define VISC_GRID
# define SPONGE
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# undef TS_MPDATA
#endif


#ifdef SOLVE3D
# define TS_DIF2
# define MIX_GEO_TS
# define DIFF_GRID
#endif


/* vertical mixing */

#ifdef SOLVE3D
# define SOLAR_SOURCE

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
# define SRELAXATION
# undef QCORRECTION
#endif

#ifdef SOLVE3D
# undef TCLIMATOLOGY
# undef TCLM_NUDGING
#endif

/* point sources (rivers, line sources) */

/* Using Runoff instead now */
#undef ANA_PSOURCE
#ifdef SOLVE3D
# define RUNOFF
#endif

/* tides */

#undef LTIDES
#ifdef LTIDES
# undef FILTERED
# define SSH_TIDES
# define UV_TIDES
# ifdef SOLVE3D
#  define ADD_FSOBC
#  define ADD_M2OBC
# endif
# undef RAMP_TIDES
# define TIDES_ASTRO
# undef POT_TIDES

# undef UV_LDRAG
# define UV_DRAG_GRID
# define ANA_DRAG
# define LIMIT_BSTRESS
# define UV_QDRAG
#else
# undef UV_QDRAG
# define UV_LDRAG
#endif

#define RADIATION_2D

/* roms quirks */

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# undef ANA_SMFLUX
# define BULK_FLUXES2D
# define CCSM_FLUXES2D
# define ANA_INITIAL
# define ANA_FSOBC
# define ANA_M2OBC
#endif
