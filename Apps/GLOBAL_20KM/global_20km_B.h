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
**  Options for dynamical downscaling of IPCC AR4 GISS AOM results
**  to the 10 km Atlantic-Arctic grid
*/

#if defined GLOBAL_20KM
!#define INLINE_2DIO        /* define if processing 3D IO level by level */
!#define NO_HIS              /* define if no history files to be created or written */
#define RST_SINGLE
#undef  PERFECT_RESTART
#define ERA40               /* 20C3M 1981-2000 but with ERA40 */
#undef  STATIONS
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#undef  VISC_3DCOEF
#undef  UV_SMAGORINSKY
!#define UV_U3ADV_SPLIT
#define MIX_S_UV
!#define MIX_GEO_UV
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
!#define TS_DIF2
#undef  TS_DIF4
#define MIX_GEO_TS
!#define DIFF_3DCOEF
!#define TS_SMAGORINSKY
!#define TS_U3ADV_SPLIT
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#define MASKING
#define AVERAGES
#define UV_QDRAG        /* turn ON or OFF quadratic bottom friction */
#undef  MY25_MIXING
#undef  GLS_MIXING
#ifdef GLS_MIXING
#  define  N2S2_HORAVG
#endif
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
# define LMD_DDMIX
#endif
#define ANA_BSFLUX
#define ANA_BTFLUX
#undef  ANA_GRID
#undef  ANA_INITIAL
#undef  ANA_MASK
#undef  ANA_MEANRHO
#define ANA_SSFLUX
#define ANA_STFLUX
#undef  ANA_SMFLUX
#undef  ANA_SRFLUX
#undef  ANA_VMIX
#undef  ANA_CLOUD
#undef  ANA_AIRT
#undef  ANA_DEWT
#undef  ANA_SLP
#undef  ANA_WINDS

#define EASTERN_WALL
#define WESTERN_WALL
#define NS_PERIODIC

#define SCORRECTION

#define ICE_MODEL
# ifdef ICE_MODEL
#  define ICE_THERMO
#    define ICE_MK
#    undef ICE_ALB_EC92
#  define ICE_MOMENTUM
#  define ICE_BULK_FLUXES
#    undef  ICE_MOM_BULK
#    define ICE_EVP
#  define ICE_ADVECT
#    define ICE_SMOLAR
#    define ICE_UPWIND
# endif

#define SPECIFIC_HUMIDITY

#define BULK_FLUXES
#ifdef  BULK_FLUXES
# ifdef ICE_MODEL
#  define ICE_BULK_FLUXES
# endif
# define EMINUSP
# define ALBEDO_CLOUD
# ifdef ERA40
#  define CLOUDS
#  define LONGWAVE
#  define COOL_SKIN
#  define ANA_SRFLUX
#  define SHORTWAVE
# else
#  undef  CLOUDS
#  define LONGWAVE_OUT
#  define COOL_SKIN
#  undef  ANA_SRFLUX
#  define SHORTWAVE
# endif
#endif

#undef  NCEP_FLUXES

#define RUNOFF

#define SOLAR_SOURCE
#undef  SLP_GRAD

#undef   SSSFLX

#endif
