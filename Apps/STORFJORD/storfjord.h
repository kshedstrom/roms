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

#if defined STORFJORD
!#define INLINE_2DIO        /* define if processing 3D IO level by level */
!#define NO_HIS              /* define if no history files to be created or written */
!#define RST_SINGLE
#define PERFECT_RESTART
#define NO_WRITE_GRID
#define ERA_INTERIM
#define STATIONS
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#define UV_VIS2
#undef  UV_VIS4
#define VISC_3DCOEF
#define UV_SMAGORINSKY
#undef  UV_U3ADV_SPLIT
#define MIX_S_UV
#undef  MIX_GEO_UV
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  WJ_GRADP
#define DJ_GRADPS
!#define  DIFF_GRID
!#define TS_DIF2
!#define  TS_DIF4
!#define MIX_GEO_TS
!#define MIX_S_TS
!#define DIFF_3DCOEF
!#define TS_SMAGORINSKY
#undef  TS_U3ADV_SPLIT
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
!#define TS_C4HADVECTION
!#define TS_C4VADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#undef  TS_MPDATA
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#define MASKING
#define AVERAGES
!#define AVERAGES_DETIDE
#undef  UV_QDRAG        /* turn ON or OFF quadratic bottom friction */
#define UV_LDRAG        /* turn ON or OFF lineat bottom friction */
#define UV_DRAG_GRID
#define DRAG_LIMITER
#undef  MY25_MIXING
#define GLS_MIXING
#ifdef GLS_MIXING
#  define  N2S2_HORAVG
#endif
#undef  LMD_MIXING
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

#define TCLIMATOLOGY
#define M3CLIMATOLOGY
#define TCLM_NUDGING
#define M3CLM_NUDGING

#undef  SCORRECTION


#define SSH_TIDES       /* turn on computation of tidal elevation */
#define UV_TIDES        /* turn on computation of tidal currents */
#define ADD_FSOBC       /* Add tidal elevation to processed OBC data */
#define ADD_M2OBC       /* Add tidal currents  to processed OBC data */
#define RAMP_TIDES      /* Spin up tidal forcing */
#define TIDES_ASTRO     /* apply nodal corrections */
#undef  POT_TIDES        /* turn on computation of tidal potential */

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

#define CCSM_FLUXES
#define BULK_FLUXES
#ifdef  BULK_FLUXES
# ifdef ICE_MODEL
#  define ICE_BULK_FLUXES
# endif
# define EMINUSP
# define ALBEDO_CLOUD
# ifdef ERA_INTERIM
#  define CLOUDS
#  define LONGWAVE_OUT
#  define COOL_SKIN
#  undef  ANA_SRFLUX
#  define SHORTWAVE
# elif defined ERA40
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

#undef  RUNOFF

#define SOLAR_SOURCE
#undef  SLP_GRAD

#undef   SSSFLX

#endif
