#if defined SMALL_AREA
/*
  Options for Atlantic + Arctic model
*/
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef LMD_BKPP
# define LMD_NONLOCAL
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
#undef  TCLM_NUDGING
#undef  TCLIMATOLOGY
#define EASTERN_WALL
#define NORTHERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#undef  EAST_FSGRADIENT
#undef  EAST_M2GRADIENT
#undef  EAST_M3RADIATION
#undef  EAST_TRADIATION
!
#define NCEP_FLUXES
#define SOLAR_SOURCE
!
#define ICE_MODEL
# define ICE_THERMO
#   undef ICE_SMOOTH
# define ICE_MOMENTUM
#   define ICE_MOM_BULK
#   define ICE_EVP
# define ICE_ADVECT
#   define ICE_SMOLAR
#     define ICE_UPWIND
!
# elif defined LARGE_AREA
/*
  Options for Atlantic + Arctic model
*/
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef LMD_BKPP
# define LMD_NONLOCAL
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
!
#define OBC_VOLCONS
#define ANA_M2OBC
#define WEST_M2CLAMPED
#define WEST_M3CLAMPED
#define EAST_M2CLAMPED
#define EAST_M3CLAMPED
#define ANA_M3OBC
#define TCLM_NUDGING
#define TCLIMATOLOGY
#undef  EASTERN_WALL
#undef  WESTERN_WALL
#define NORTHERN_WALL
#define SOUTHERN_WALL
#undef EAST_FSGRADIENT
#undef EAST_M2GRADIENT
#undef EAST_M3RADIATION
#undef EAST_TRADIATION
!
#define NCEP_FLUXES
#define SOLAR_SOURCE
!
#define ICE_MODEL
# define ICE_THERMO
#   undef ICE_ALB_EC92
# define ICE_MOMENTUM
#   undef  ICE_MOM_BULK
#   define ICE_EVP
# define ICE_ADVECT
#   define ICE_SMOLAR
#     define ICE_UPWIND
!
# elif defined LARGE_AREA_0_8
/*
  Options for Atlantic + Arctic model
*/
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef LMD_BKPP
# define LMD_NONLOCAL
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
!
#define OBC_VOLCONS
#define ANA_M2OBC
#define WEST_M2CLAMPED
#define WEST_M3CLAMPED
#define EAST_M2CLAMPED
#define EAST_M3CLAMPED
#define ANA_M3OBC
#define TCLM_NUDGING
#define TCLIMATOLOGY
#undef  EASTERN_WALL
#undef  WESTERN_WALL
#define NORTHERN_WALL
#define SOUTHERN_WALL
#undef EAST_FSGRADIENT
#undef EAST_M2GRADIENT
#undef EAST_M3RADIATION
#undef EAST_TRADIATION
!
#define NCEP_FLUXES
#define SOLAR_SOURCE
!
#define ICE_MODEL
# define ICE_THERMO
#   undef ICE_ALB_EC92
# define ICE_MOMENTUM
#   undef  ICE_MOM_BULK
#   define ICE_EVP
# define ICE_ADVECT
#   define ICE_SMOLAR
#     define ICE_UPWIND
!
# elif defined BARENTS_0_2
/*
  Options for Barents Sea model
*/
#define INLINE_2DIO        /* define if processing 3D IO level by level */
#define RST_SINGLE
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#define UV_SADVECTION
#define UV_QDRAG        /* turn ON or OFF quadratic bottom friction */
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#undef  MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#define GLS_MIXING
#ifdef GLS_MIXING
#  define  N2S2_HORAVG
#endif
#undef LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
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
!
#undef  OBC_VOLCONS
#undef  ANA_M2OBC
#undef  ANA_M3OBC
!
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_M3FRS
#define WEST_TFRS
#define WEST_AICLAMPED
#define WEST_HICLAMPED
#define WEST_HSNCLAMPED
#define WEST_TICLAMPED
#define WEST_SFWATCLAMPED
#define WEST_SIG11CLAMPED
#define WEST_SIG22CLAMPED
#define WEST_SIG12CLAMPED
#define WEST_MIGRADIENT
!
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3FRS
#define EAST_TFRS
#define EAST_AICLAMPED
#define EAST_HICLAMPED
#define EAST_HSNCLAMPED
#define EAST_TICLAMPED
#define EAST_SFWATCLAMPED
#define EAST_SIG11CLAMPED
#define EAST_SIG22CLAMPED
#define EAST_SIG12CLAMPED
#define EAST_MIGRADIENT
!
#define NORTH_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_M3FRS
#define NORTH_TFRS
#define NORTH_AICLAMPED
#define NORTH_HICLAMPED
#define NORTH_HSNCLAMPED
#define NORTH_TICLAMPED
#define NORTH_SFWATCLAMPED
#define NORTH_SIG11CLAMPED
#define NORTH_SIG22CLAMPED
#define NORTH_SIG12CLAMPED
#define NORTH_MIGRADIENT
!
#define SOUTHERN_WALL
!
!
#define  TCLIMATOLOGY
#define  M3CLIMATOLOGY
!
/* Forcing */
#define  SSH_TIDES       /* turn on computation of tidal elevation */
#define  UV_TIDES        /* turn on computation of tidal currents */
#define  ADD_FSOBC       /* Add tidal elevation to processed OBC data */
#define  ADD_M2OBC       /* Add tidal currents  to processed OBC data */
#define  TIDES_NOSPIN    /* Do not spin up tidal forcing */
!
#define ICE_MODEL
# ifdef ICE_MODEL
#  define ICE_THERMO
#    define ICE_MK
#    define ICE_ALB_EC92
#  define ICE_MOMENTUM
#    undef  ICE_MOM_BULK
#    define ICE_EVP
#  define ICE_ADVECT
#    define ICE_SMOLAR
#    define ICE_UPWIND
# endif
!
#define  BULK_FLUXES
#ifdef  BULK_FLUXES
# ifdef ICE_MODEL
#  define ICE_BULK_FLUXES
# endif
# define EMINUSP
# define ANA_SRFLUX
# define ALBEDO
# define CLOUDS
# define LONGWAVE
# define SHORTWAVE
#endif
!
#undef NCEP_FLUXES
!
#define RUNOFF
!
#define SOLAR_SOURCE
!
/* Rivers */
#undef UV_PSOURCE
#undef TS_PSOURCE
!
# elif defined BARENTS_1D
/*
  Options for idealized 1-D ice-ocean system
*/
#undef  ANA_GRID
#undef  ANA_INITIAL
#undef  DIAGNOSTIC
#undef  CURVGRID
#define UV_ADV
#define UV_QDRAG        /* turn ON or OFF quadratic bottom friction */
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#undef  UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#undef  TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
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
# undef LMD_BKPP
# define LMD_NONLOCAL
#endif
#define ANA_BSFLUX
#define ANA_BTFLUX
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
#undef  TCLM_NUDGING
#undef  TCLIMATOLOGY
#undef  EASTERN_WALL
#undef  NORTHERN_WALL
#undef  WESTERN_WALL
#undef  SOUTHERN_WALL
#undef  EAST_FSGRADIENT
#undef  EAST_M2GRADIENT
#undef  EAST_M3RADIATION
#undef  EAST_TRADIATION
#define NS_PERIODIC
#define EW_PERIODIC
!
#define ICE_MODEL
# undef ANA_ICE
# define ICE_THERMO
# ifdef ICE_THERMO
#   define ICE_MK
#   define ICE_ALB_EC92
#   undef ICE_SMOOTH
# endif
# define ICE_MOMENTUM
#   undef  ICE_MOM_BULK
#   define ICE_EVP
# define ICE_ADVECT
#   define ICE_SMOLAR
#     define ICE_UPWIND
!
#define BULK_FLUXES
#ifdef  BULK_FLUXES
# ifdef ICE_MODEL
#  define ICE_BULK_FLUXES
# endif
# define EMINUSP
# define ANA_SRFLUX
# define ALBEDO
# define CLOUDS
# define LONGWAVE
# define SHORTWAVE
#endif
!
#undef  NCEP_FLUXES
# undef ANA_NCEP
!
#define SOLAR_SOURCE
!
!
# elif defined CONMAN
/*
  Options for CONMAN 4-km grid model
*/
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#undef  LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
#endif
#define GLS_MIXING
#ifdef GLS
# define N2S2_HORAVG
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
!
#undef  ANA_M2OBC
#undef  ANA_M3OBC
!
!
#undef  OBC_VOLCONS
!
#undef  WESTERN_WALL
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_M3FRS
#define WEST_M3CLAMPED
#define WEST_TFRS
#define WEST_TCLAMPED
#undef WEST_TRADIATION
#undef WEST_TNUDGING
!
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3FRS
#define EAST_M3CLAMPED
#define EAST_TFRS
#define EAST_TCLAMPED
#undef   EAST_TRADIATION
#undef   EAST_TNUDGING
!
#define  SOUTHERN_WALL
!
#define NORTH_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_M3FRS
#define NORTH_M3CLAMPED
#define NORTH_TFRS
#define NORTH_TCLAMPED
#undef   NORTH_TRADIATION
#undef   NORTH_TNUDGING
!
!
#define  TCLIMATOLOGY
#define  M3CLIMATOLOGY
!
/* Forcing */
#define  SSH_TIDES       /* turn on computation of tidal elevation */
#define  UV_TIDES        /* turn on computation of tidal currents */
#define  ADD_FSOBC       /* Add tidal elevation to processed OBC data */
#define  ADD_M2OBC       /* Add tidal currents  to processed OBC data */
!
#define NCEP_FLUXES
#define SOLAR_SOURCE
!
/* Rivers */
!
#undef   ANA_PSOURCE
#undef UV_PSOURCE
#undef TS_PSOURCE
!
/* write station data */
#define STATIONS
!
#undef ICE_MODEL
!
!
# elif defined MOZAMBIQUE_9KM
/*
  Options for Mozambique regional 9-km model
*/
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
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
!
#define UV_PSOURCE
#define TS_PSOURCE
#define ANA_PSOURCE
!
#define OBC_VOLCONS
#undef  ANA_M2OBC
#undef  ANA_M3OBC
!
!
#define OBC_VOLCONS
#undef  WESTERN_WALL
#undef   WEST_FSGRADIENT
#undef   WEST_M2GRADIENT
#define WEST_M3GRADIENT
#define WEST_TRADIATION
#undef   SOUTHERN_WALL
#undef   SOUTH_FSGRADIENT
#undef   SOUTH_M2GRADIENT
#define SOUTH_M3GRADIENT
#define SOUTH_TRADIATION
#define NORTHERN_WALL
#define  EASTERN_WALL
#undef   EAST_FSGRADIENT
#undef   EAST_M2GRADIENT
#undef   EAST_M3GRADIENT
#undef   EAST_TRADIATION
!
#undef  TCLM_NUDGING
#define TCLIMATOLOGY
!
!
!
#define NCEP_FLUXES
#define SOLAR_SOURCE
!
#undef ICE_MODEL
!
# elif defined WAP_0_1
/*
  Options for Western Antarctic Peninsula model
*/
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef  LMD_BKPP
# define LMD_NONLOCAL
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
!
#define WEST_VOLCONS
#define EAST_VOLCONS
#define NORTH_VOLCONS
#define OBC_VOLCONS
#undef  ANA_M2OBC
#undef  ANA_M3OBC
!
#undef  WESTERN_WALL
#define WEST_M2RADIATION
#define WEST_M3RADIATION
#define WEST_TRADIATION
#define WEST_AIGRADIENT
#define WEST_HIGRADIENT
#define WEST_HSNGRADIENT
#define WEST_TIGRADIENT
#define WEST_SFWATGRADIENT
#define WEST_SIG11GRADIENT
#define WEST_SIG22GRADIENT
#define WEST_SIG12GRADIENT
#define WEST_MIGRADIENT
!
#define  EASTERN_WALL
#ifndef EASTERN_WALL
#define EAST_M2GRADIENT
#define EAST_M3GRADIENT
#define EAST_TRADIATION
#define EAST_AIGRADIENT
#define EAST_HIGRADIENT
#define EAST_HSNGRADIENT
#define EAST_TIGRADIENT
#define EAST_SFWATGRADIENT
#define EAST_SIG11GRADIENT
#define EAST_SIG22GRADIENT
#define EAST_SIG12GRADIENT
#define EAST_MIGRADIENT
#endif
!
#undef  NORTHERN_WALL
#define NORTH_M2RADIATION
#define NORTH_M3RADIATION
#define NORTH_TRADIATION
#define NORTH_AIGRADIENT
#define NORTH_HIGRADIENT
#define NORTH_HSNGRADIENT
#define NORTH_TIGRADIENT
#define NORTH_SFWATGRADIENT
#define NORTH_SIG11GRADIENT
#define NORTH_SIG22GRADIENT
#define NORTH_SIG12GRADIENT
#define NORTH_MIGRADIENT
!
#define SOUTHERN_WALL
!
!
#undef  TCLM_NUDGING
#define TCLIMATOLOGY
!
!
!
#define NCEP_FLUXES
#define SOLAR_SOURCE
!
#define ICE_MODEL
#ifdef ICE_MODEL
# define ICE_THERMO
#   undef ICE_ALB_EC92
# define ICE_MOMENTUM
#   undef  ICE_MOM_BULK
#   define ICE_EVP
# define ICE_ADVECT
#   define ICE_SMOLAR
#     define ICE_UPWIND
#endif
!
!
# elif defined WAP_0_1_1
/*
  Options for Western Antarctic Peninsula model larger domain
*/
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#define UV_VIS2
#undef  UV_VIS4
#define UV_SADVECTION
#define VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef  LMD_BKPP
# define LMD_NONLOCAL
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
!
#define RADIATION_2D
!
#define WEST_VOLCONS
#define EAST_VOLCONS
#define NORTH_VOLCONS
#define OBC_VOLCONS
#undef  ANA_M2OBC
#undef  ANA_M3OBC
!
#undef  WESTERN_WALL
#define WEST_FSRADIATION
#define WEST_M2RADIATION
#define WEST_M3RADIATION
#define WEST_TRADIATION
#define WEST_AIGRADIENT
#define WEST_HIGRADIENT
#define WEST_HSNGRADIENT
#define WEST_TIGRADIENT
#define WEST_SFWATGRADIENT
#define WEST_SIG11GRADIENT
#define WEST_SIG22GRADIENT
#define WEST_SIG12GRADIENT
#define WEST_MIGRADIENT
!
#undef  EASTERN_WALL
#ifndef EASTERN_WALL
#define EAST_FSRADIATION
#define EAST_M2RADIATION
#define EAST_M3RADIATION
#define EAST_TRADIATION
#define EAST_AIGRADIENT
#define EAST_HIGRADIENT
#define EAST_HSNGRADIENT
#define EAST_TIGRADIENT
#define EAST_SFWATGRADIENT
#define EAST_SIG11GRADIENT
#define EAST_SIG22GRADIENT
#define EAST_SIG12GRADIENT
#define EAST_MIGRADIENT
#endif
!
#undef  NORTHERN_WALL
#define NORTH_FSRADIATION
#define NORTH_M2RADIATION
#define NORTH_M3RADIATION
#define NORTH_TRADIATION
#define NORTH_AIGRADIENT
#define NORTH_HIGRADIENT
#define NORTH_HSNGRADIENT
#define NORTH_TIGRADIENT
#define NORTH_SFWATGRADIENT
#define NORTH_SIG11GRADIENT
#define NORTH_SIG22GRADIENT
#define NORTH_SIG12GRADIENT
#define NORTH_MIGRADIENT
!
#define SOUTHERN_WALL
!
!
#undef  TCLM_NUDGING
#define TCLIMATOLOGY
!
!
!
#define NCEP_FLUXES
#define SOLAR_SOURCE
!
#define ICE_MODEL
#ifdef ICE_MODEL
# define ICE_THERMO
#   undef ICE_ALB_EC92
# define ICE_MOMENTUM
#   undef  ICE_MOM_BULK
#   define ICE_EVP
# define ICE_ADVECT
#   define ICE_SMOLAR
#     define ICE_UPWIND
#endif
# elif defined WAP_0_1_2
/*
  Options for Western Antarctic Peninsula model larger domain
*/
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#define UV_VIS2
#undef  UV_VIS4
#define UV_SADVECTION
#define VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef  LMD_BKPP
# define LMD_NONLOCAL
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
!
#define RADIATION_2D
!
#define WEST_VOLCONS
#define EAST_VOLCONS
#define NORTH_VOLCONS
#define OBC_VOLCONS
#undef  ANA_M2OBC
#undef  ANA_M3OBC
!
#undef  WESTERN_WALL
#define WEST_FSRADIATION
#define WEST_M2RADIATION
#define WEST_M3RADIATION
#define WEST_TRADIATION
#define WEST_AIGRADIENT
#define WEST_HIGRADIENT
#define WEST_HSNGRADIENT
#define WEST_TIGRADIENT
#define WEST_SFWATGRADIENT
#define WEST_SIG11GRADIENT
#define WEST_SIG22GRADIENT
#define WEST_SIG12GRADIENT
#define WEST_MIGRADIENT
!
#undef  EASTERN_WALL
#ifndef EASTERN_WALL
#define EAST_FSRADIATION
#define EAST_M2RADIATION
#define EAST_M3RADIATION
#define EAST_TRADIATION
#define EAST_AIGRADIENT
#define EAST_HIGRADIENT
#define EAST_HSNGRADIENT
#define EAST_TIGRADIENT
#define EAST_SFWATGRADIENT
#define EAST_SIG11GRADIENT
#define EAST_SIG22GRADIENT
#define EAST_SIG12GRADIENT
#define EAST_MIGRADIENT
#endif
!
#undef  NORTHERN_WALL
#define NORTH_FSRADIATION
#define NORTH_M2RADIATION
#define NORTH_M3RADIATION
#define NORTH_TRADIATION
#define NORTH_AIGRADIENT
#define NORTH_HIGRADIENT
#define NORTH_HSNGRADIENT
#define NORTH_TIGRADIENT
#define NORTH_SFWATGRADIENT
#define NORTH_SIG11GRADIENT
#define NORTH_SIG22GRADIENT
#define NORTH_SIG12GRADIENT
#define NORTH_MIGRADIENT
!
#define SOUTHERN_WALL
!
!
#undef  TCLM_NUDGING
#define TCLIMATOLOGY
!
!
!
#define NCEP_FLUXES
#define SOLAR_SOURCE
!
#define ICESHELF
!
#define ICE_MODEL
#ifdef ICE_MODEL
# define ICE_THERMO
#   undef ICE_ALB_EC92
# define ICE_MOMENTUM
#   undef  ICE_MOM_BULK
#   define ICE_EVP
# define ICE_ADVECT
#   define ICE_SMOLAR
#     define ICE_UPWIND
#endif
!
# elif defined ACC_1
/*
  Options for ACC experiment
*/
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#define UV_QDRAG
#undef  UV_VIS2
#undef  MIX_S_UV
#undef  SMAGORINSKY
#undef  UV_VIS4
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef  LMD_BKPP
# define LMD_NONLOCAL
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
!
!
#define EW_PERIODIC
!
#define NORTHERN_WALL
!
#define SOUTHERN_WALL
!
!
#undef  TCLM_NUDGING
#define TCLIMATOLOGY
!
!
!
#define NCEP_FLUXES
#define SOLAR_SOURCE
!
#define ICESHELF
!
#define ICE_MODEL
#ifdef ICE_MODEL
# define ICE_THERMO
#   define ICE_MK
#   undef ICE_ALB_EC92
# define ICE_MOMENTUM
#   undef  ICE_MOM_BULK
#   define ICE_EVP
# define ICE_ADVECT
#   define ICE_SMOLAR
#     define ICE_UPWIND
#endif
!
# elif defined GAK1D
#define UV_ADV
#define UV_VIS2
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define TS_DIF2
#define SALINITY
#define NONLIN_EOS
#define SOLAR_SOURCE
#define DJ_GRADPS
#define BULK_FLUXES
#define LONGWAVE
#define ALBEDO
#define SOLVE3D
#define PROFILE
#define AVERAGES
#define AVERAGES_AKS
#define AVERAGES_FLUXES
#define SPHERICAL
#define SPLINES
#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_SRFLUX
#define LMD_MIXING
# ifdef LMD_MIXING
#define LMD_RIMIX
#define LMD_CONVEC
#undef  LMD_DDMIX
#define LMD_SKPP
#define LMD_BKPP
#define LMD_NONLOCAL
# endif
#define EW_PERIODIC
#define NS_PERIODIC
#define BIO_GOANPZ
!
  /*
**  Biological model options.
*/
!
#undef  NPZD1              /* Craig Lewiss 4 box model */
#define BIO_GOANPZ         /* Sarah Hinckleys 10 box model */
#undef  PASSIVE_TRACERS    /* add 5 tracer boxes that are passive */
!
#if defined NPZD1 || defined BIO_GOANPZ || defined PASSIVE_TRACERS
# define BIOFLUX           /* sum Nitrogen fluxes between boxes */
# define ANA_BIOLOGY       /* analytical biology initial conditions */
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
!
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
#endif
#ifdef BIO_GOANPZ
# define OFFLINE_BIOLOGY   /* define if offline simulation of bio tracers */
# define DIAPAUSE          /* Enable Neocalanus seasonal vertical migration */
# define IRON_LIMIT              /* Add iron as passive 13th tracer */
# if defined IRON_LIMIT
#  if !defined OFFLINE_BIOLOGY
#  define  ANA_TCLIMA      /* analytical tracers climatology for iron */
#  define  TCLIMATOLOGY    /* Processing of tracer climatology for iron */
#  endif
#  define  TCLM_NUDGING    /* Nudging of tracer climatology for iron */
# endif
#endif
!
!
# elif defined ICE_BASIN
/*
  Options for idealized ice-covered basin
*/
#define ANA_GRID
#define ANA_INITIAL
#undef  DIAGNOSTIC
#undef  CURVGRID
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#undef  UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#undef TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef LMD_BKPP
# define LMD_NONLOCAL
#endif
#define ANA_BSFLUX
#define ANA_BTFLUX
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
#undef  TCLM_NUDGING
#undef  TCLIMATOLOGY
#define EASTERN_WALL
#define NORTHERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#undef  EAST_FSGRADIENT
#undef  EAST_M2GRADIENT
#undef  EAST_M3RADIATION
#undef  EAST_TRADIATION
!
#define NCEP_FLUXES
# define ANA_NCEP
#define SOLAR_SOURCE
!
#define ICE_MODEL
# define ANA_ICE
# define ICE_THERMO
#   undef ICE_SMOOTH
# define ICE_MOMENTUM
#   undef  ICE_MOM_BULK
#   define ICE_EVP
# define ICE_ADVECT
#   define ICE_SMOLAR
#     define ICE_UPWIND
!
# elif defined ICE_OCEAN_1D
/*
  Options for idealized 1-D ice-ocean system
*/
#define ANA_GRID
#define ANA_INITIAL
#undef  DIAGNOSTIC
#undef  CURVGRID
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#undef  UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#undef TS_SVADVECTION
#define MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef LMD_BKPP
# define LMD_NONLOCAL
#endif
#define ANA_BSFLUX
#define ANA_BTFLUX
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
#undef  TCLM_NUDGING
#undef  TCLIMATOLOGY
#undef  EASTERN_WALL
#undef  NORTHERN_WALL
#undef  WESTERN_WALL
#undef  SOUTHERN_WALL
#undef  EAST_FSGRADIENT
#undef  EAST_M2GRADIENT
#undef  EAST_M3RADIATION
#undef  EAST_TRADIATION
#define NS_PERIODIC
#define EW_PERIODIC
!
#define NCEP_FLUXES
# define ANA_NCEP
#define SOLAR_SOURCE
!
#define ICE_MODEL
# define ANA_ICE
# define ICE_THERMO
#   undef ICE_SMOOTH
# define ICE_MOMENTUM
#   undef  ICE_MOM_BULK
#   define ICE_EVP
# define ICE_ADVECT
#   define ICE_SMOLAR
#     define ICE_UPWIND
!
# elif defined GLOBAL_20KM
/*
  Options for Global Ice-Ocean Model
*/
#define INLINE_2DIO        /* define if processing 3D IO level by
level */
#define RST_SINGLE
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#define UV_VIS2
#define SMAGORINSKY
#undef  UV_VIS4
#define UV_QDRAG
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#define MIX_S_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#undef  MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#define STATIONS
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
# undef  LMD_BKPP
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
!
#define NS_PERIODIC
!
#undef  TCLM_NUDGING
#undef  TCLIMATOLOGY
#define SCORRECTION
!
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
#    undef  ICE_UPWIND
# endif
!
#define  BULK_FLUXES
#ifdef  BULK_FLUXES
# ifdef ICE_MODEL
#  define ICE_BULK_FLUXES
# endif
# define EMINUSP
# define ALBEDO
# undef  CLOUDS
# define LONGWAVE_OUT
# define COOL_SKIN
# undef  ANA_SRFLUX
# define SHORTWAVE
#endif
!
#undef  NCEP_FLUXES
!
#define RUNOFF
!
#define SOLAR_SOURCE
#define SLP_GRAD
!
#define POT_TIDES       /* turn on computation of tidal potential
*/
#define TIDES_ASTRO     /* apply nodal corrections */
#define TIDES_NOSPIN
!
#define  WRT_SSSFLX
#undef   SSSFLX
!
!
# elif defined GLOBAL_TEST
/*
  Options for Global Ice-Ocean Model
*/
#define INLINE_2DIO        /* define if processing 3D IO level by
level */
#define RST_SINGLE
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#define UV_VIS2
#undef  UV_VIS4
#define UV_QDRAG
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#define MIX_S_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#undef  MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#undef  BODYFORCE
#define SPLINES
#define MASKING
#define AVERAGES
#define STATIONS
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
!
#define NS_PERIODIC
!
#undef   TCLM_NUDGING
#undef   TCLIMATOLOGY
#define  SCORRECTION
!
#define ICE_MODEL
# ifdef ICE_MODEL
#  define ICE_THERMO
#    define ICE_MK
#    define ICE_ALB_EC92
#  define ICE_MOMENTUM
#    undef  ICE_MOM_BULK
#    define ICE_EVP
#  define ICE_ADVECT
#    define ICE_SMOLAR
#    define ICE_UPWIND
# endif
!
#define BULK_FLUXES
#ifdef  BULK_FLUXES
# ifdef ICE_MODEL
#  define ICE_BULK_FLUXES
# endif
# define EMINUSP
# define ALBEDO
# undef  CLOUDS
# undef  LONGWAVE
# define LONGWAVE_OUT
# define SHORTWAVE
# undef  ANA_SRFLUX
#endif
!
#undef  NCEP_FLUXES
!
#define RUNOFF
!
#define SOLAR_SOURCE
#define SLP_GRAD
!
#define POT_TIDES       /* turn on computation of tidal potential
*/
#define TIDES_ASTRO     /* apply nodal corrections */
#undef  TIDES_NOSPIN
!
#define  WRT_SSSFLX
#undef   SSSFLX

# elif defined MARECO
/*
  Options for Mid-Atlantic Ridge model
*/
#define INLINE_2DIO        /* define if processing 3D IO level by level */
!
#define RST_SINGLE
#undef  DIAGNOSTIC
#define CURVGRID
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#define UV_SADVECTION
#undef  VISC_GRID
#define NONLIN_EOS
#undef  MIX_GEO_UV
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  DIFF_GRID
#undef  TS_DIF2
#undef  TS_DIF4
#define TS_U3HADVECTION
#undef  TS_A4HADVECTION
#undef  TS_C4HADVECTION
#undef  TS_A4VADVECTION
#define TS_SVADVECTION
#define MIX_GEO_TS
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
# undef  LMD_BKPP
# define LMD_NONLOCAL
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
!
#undef  OBC_VOLCONS
#undef  ANA_M2OBC
#undef  ANA_M3OBC
!
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_M3FRS
#define WEST_TFRS
#define WEST_AICLAMPED
#define WEST_HICLAMPED
#define WEST_HSNCLAMPED
#define WEST_TICLAMPED
#define WEST_SFWATCLAMPED
#define WEST_SIG11CLAMPED
#define WEST_SIG22CLAMPED
#define WEST_SIG12CLAMPED
#define WEST_MIGRADIENT
!
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3FRS
#define EAST_TFRS
#define EAST_AICLAMPED
#define EAST_HICLAMPED
#define EAST_HSNCLAMPED
#define EAST_TICLAMPED
#define EAST_SFWATCLAMPED
#define EAST_SIG11CLAMPED
#define EAST_SIG22CLAMPED
#define EAST_SIG12CLAMPED
#define EAST_MIGRADIENT
!
#define NORTH_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_M3FRS
#define NORTH_TFRS
#define NORTH_AICLAMPED
#define NORTH_HICLAMPED
#define NORTH_HSNCLAMPED
#define NORTH_TICLAMPED
#define NORTH_SFWATCLAMPED
#define NORTH_SIG11CLAMPED
#define NORTH_SIG22CLAMPED
#define NORTH_SIG12CLAMPED
#define NORTH_MIGRADIENT
!
#define SOUTH_FSCHAPMAN
#define SOUTH_M2FLATHER
#define SOUTH_M3FRS
#define SOUTH_TFRS
#define SOUTH_AICLAMPED
#define SOUTH_HICLAMPED
#define SOUTH_HSNCLAMPED
#define SOUTH_TICLAMPED
#define SOUTH_SFWATCLAMPED
#define SOUTH_SIG11CLAMPED
#define SOUTH_SIG22CLAMPED
#define SOUTH_SIG12CLAMPED
#define SOUTH_MIGRADIENT
!
#define  TCLIMATOLOGY
#define  M3CLIMATOLOGY
!
/* Forcing */
#undef  SSH_TIDES       /* turn on computation of tidal elevation
*/
#undef  UV_TIDES        /* turn on computation of tidal currents
*/
#undef  ADD_FSOBC       /* Add tidal elevation to processed OBC
data */
#undef  ADD_M2OBC       /* Add tidal currents  to processed OBC
data */
#undef  TIDES_NOSPIN    /* Do not spin up tidal forcing */
!
#define NCEP_FLUXES
#define SOLAR_SOURCE
#undef  RUNOFF
#define SCORRECTION
!
#undef   BULK_FLUXES
#undef   ICE_BULK_FLUXES
#undef   SPECIFIC_HUMIDITY
#undef   LONGWAVE
#undef   SHORTWAVE
#undef   ANA_SRFLUX
#undef   CLOUDS
#undef   EVAPORATION
#undef   PRECIPITATION
#undef   SOLAR_SOURCE
#undef   SCORRECTION
!
/* Rivers */
#undef   UV_PSOURCE
#undef   TS_PSOURCE
!
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
#endif
