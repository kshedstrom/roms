/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2011 The ROMS/TOMS Group
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

#if defined GISS_AOM
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
#define UV_VIS2
#undef  UV_VIS4
#define VISC_3DCOEF
#define UV_SMAGORINSKY
!#define UV_U3ADV_SPLIT
#undef  MIX_S_UV
#define MIX_GEO_UV
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#define TS_DIF2
#undef  TS_DIF4
#define MIX_GEO_TS
#define DIFF_3DCOEF
#define TS_SMAGORINSKY
!#define TS_U3ADV_SPLIT
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
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

#undef  OBC_VOLCONS
#undef  ANA_M2OBC
#undef  ANA_M3OBC

#define WESTERN_WALL

#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_TRADIATION
#define EAST_TNUDGING
#define EAST_M3RADIATION
#define EAST_M3NUDGING
#define EAST_AICLAMPED
#define EAST_HICLAMPED
#define EAST_HSNCLAMPED
#define EAST_TICLAMPED
#define EAST_SFWATCLAMPED
#define EAST_SIG11CLAMPED
#define EAST_SIG22CLAMPED
#define EAST_SIG12CLAMPED
#define EAST_MIGRADIENT

#define SOUTH_FSCHAPMAN
#define SOUTH_M2FLATHER
#define SOUTH_TRADIATION
#define SOUTH_TNUDGING
#define SOUTH_M3RADIATION
#define SOUTH_M3NUDGING
#define SOUTH_AICLAMPED
#define SOUTH_HICLAMPED
#define SOUTH_HSNCLAMPED
#define SOUTH_TICLAMPED
#define SOUTH_SFWATCLAMPED
#define SOUTH_SIG11CLAMPED
#define SOUTH_SIG22CLAMPED
#define SOUTH_SIG12CLAMPED
#define SOUTH_MIGRADIENT

#define NORTH_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_TRADIATION
#define NORTH_TNUDGING
#define NORTH_M3RADIATION
#define NORTH_M3NUDGING
#define NORTH_AICLAMPED
#define NORTH_HICLAMPED
#define NORTH_HSNCLAMPED
#define NORTH_TICLAMPED
#define NORTH_SFWATCLAMPED
#define NORTH_SIG11CLAMPED
#define NORTH_SIG22CLAMPED
#define NORTH_SIG12CLAMPED
#define NORTH_MIGRADIENT

#define TCLIMATOLOGY
#define M3CLIMATOLOGY
#define TCLM_NUDGING
#define M3CLM_NUDGING
#define NUDGING_COFF
#define SCORRECTION

/* Forcing */
#define SSH_TIDES       /* turn on computation of tidal elevation */
#define UV_TIDES        /* turn on computation of tidal currents */
#define ADD_FSOBC       /* Add tidal elevation to processed OBC data */
#define ADD_M2OBC       /* Add tidal currents  to processed OBC data */
#define RAMP_TIDES      /* Spin up tidal forcing */
#define TIDES_ASTRO     /* apply nodal corrections */
#define POT_TIDES        /* turn on computation of tidal potential */

#define ICE_MODEL
# ifdef ICE_MODEL
#  define ICE_THERMO
#    define ICE_MK
#    undef ICE_ALB_EC92
#  define ICE_MOMENTUM
#    undef  ICE_MOM_BULK
#    define ICE_EVP
#  define ICE_ADVECT
#    define ICE_SMOLAR
#    define ICE_UPWIND
# endif

#define SPECIFIC_HUMIDITY

#define  BULK_FLUXES
#ifdef  BULK_FLUXES
# ifdef ICE_MODEL
#  define ICE_BULK_FLUXES
# endif
# define EMINUSP
# define ALBEDO
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
#define SLP_GRAD

#undef   SSSFLX

#endif
