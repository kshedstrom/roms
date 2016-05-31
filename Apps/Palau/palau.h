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
**  Options for Northeast Pacific (NEP6) simulation
*/

#undef NO_HIS
#undef NESTING
#undef ONE_WAY
#define HISTORY2
#define HDF5
#define DEFLATE
#undef PARALLEL_IO
#undef OFFLINE_FLOATS

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
#undef STATIONS
#undef WET_DRY

/* output stuff */

#define NO_WRITE_GRID
#undef OUT_DOUBLE
#define PERFECT_RESTART
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

#define T_PASSIVE
#ifdef T_PASSIVE
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX 
# define ANA_PASSIVE
# undef TRC_PSOURCE
# undef ANA_TRC_PSOURCE
# define AGE_MEAN
#endif

/* advection, dissipation, pressure grad, etc. */

#ifdef SOLVE3D
# define DJ_GRADPS
#endif

#define UV_ADV
#define UV_COR
#undef UV_SADVECTION
#undef UV_C4ADVECTION

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
# define WTYPE_GRID

# define LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  define LMD_BKPP
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
#  undef DIURNAL_SRFLUX
#  define EMINUSP
#  undef ANA_SRFLUX
#  undef ALBEDO_CLOUD
#  define SOLAR_SOURCE
#  undef ALBEDO_CURVE
#  define ALBEDO_FILE
#  undef LONGWAVE
# endif
#endif

/* surface and side corrections */

#ifdef SOLVE3D
# ifdef SALINITY
#  define SCORRECTION
# else
#  define ANA_INITIAL
#  define ANA_TOBC
#  define ANA_FSOBC
#  define ANA_M2OBC
#  define ANA_M3OBC
#  define ANA_SMFLUX
#  define ANA_STFLUX
# endif
# undef QCORRECTION
#endif

/* tides */

#undef LTIDES
#ifdef LTIDES
# undef FILTERED
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# undef RAMP_TIDES
# undef TIDES_ASTRO
# undef POT_TIDES

# undef UV_LDRAG
# define UV_DRAG_GRID
# define ANA_DRAG
# define LIMIT_BSTRESS
#else
# undef M2TIDE_DIFF
#endif
#define UV_QDRAG

/* point sources (rivers, line sources) */

/* Using Runoff instead now */
#ifdef SOLVE3D
# undef RUNOFF /* no runoff in domain */
# undef ANA_PSOURCE
#endif

#define RADIATION_2D

/* roms quirks */

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif
