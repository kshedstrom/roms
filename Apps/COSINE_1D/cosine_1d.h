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
#undef ANA_INITIAL
#define ANA_GRID
#define SPHERICAL

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

#define UV_ADV
#define UV_COR
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
# define BULK_FLUXES
# define CCSM_FLUXES
# if defined BULK_FLUXES
#  define LONGWAVE_OUT
#  define SOLAR_SOURCE
#  define EMINUSP
#  undef ANA_LRFLUX
#  undef ANA_SRFLUX
#  undef ANA_ALBEDO
#  define ALBEDO
#  undef ALBEDO_CLOUD
#  undef ALBEDO_CURVE  /* for water */
#  undef LONGWAVE
#  undef ANA_PAIR
#  undef ANA_TAIR
#  undef ANA_HUMIDITY
#  undef ANA_WINDS
#  undef ANA_CLOUD
#  undef ANA_RAIN
# else
#  define ANA_SRFLUX
#  define ANA_STFLUX
#  define ANA_SSFLUX
#  define ANA_SMFLUX
# endif
#endif

#define UV_QDRAG

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

/* roms quirks */

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif

#define BIO_UMAINE
#ifdef BIO_UMAINE
# define OXYGEN
# define CARBON
# define SINK_OP2
# define TALK_NONCONSERV
# define OPTIC_MANIZZA
# undef OPTIC_UMAINE
# define IRON_LIMIT
# undef IRON_RELAX
# define ANA_BPFLUX
# define ANA_SPFLUX
#endif
