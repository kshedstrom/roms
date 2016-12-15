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
**  Options for Bneguela (BEN1) simulation
*/
 
/* general */

#define CURVGRID
#define MASKING
#define NONLIN_EOS
#define SOLVE3D
#define SALINITY
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define RI_SPLINES
# undef FLOATS
# define STATIONS
#define CCSM_COUPLED
#define CORE_FORCING
#define CCSM_FLUXES
 
/* output stuff */
 
#define NO_WRITE_GRID
#undef OUT_DOUBLE
#undef RST_SINGLE
#undef PERFECT_RESTART
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
#undef UV_SADVECTION
 
#ifdef SOLVE3D
# define TS_U3HADVECTION
# define TS_C4VADVECTION
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
# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif
#endif
 
/* surface forcing */
 
#ifdef SOLVE3D
# define BULK_FLUXES
# ifdef BULK_FLUXES
#  define LONGWAVE_OUT
#  undef DIURNAL_SRFLUX
#  define EMINUSP
#  undef ANA_SRFLUX
#  undef ALBEDO_CURVE
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
 
/* Boundary conditions...careful with grid orientation */
 
#define RADIATION_2D
#define GLOBAL_PERIODIC
 
/* roms quirks */
 
#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif
