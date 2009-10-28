#if defined LUMP

/*
**  Options for 2-d Eddy Example:
*/

#undef STATIONS
#define UV_ADV
#define UV_COR
#define UV_VIS2
#undef FLOATS
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define EW_PERIODIC
#define NS_PERIODIC

# elif defined MEDDY

/*
**  Options for Mediterranean Eddy Example:
*/

#undef STATIONS
#define UV_ADV
#define UV_COR
#undef TS_A4HADVECTION
#undef TS_A4VADVECTION
#define  TS_U3HADVECTION
#define DJ_GRADPS
#define UV_VIS2
#define TS_DIF2
#define MIX_S_UV
#define MIX_S_TS
#define SOLVE3D
#define FLOATS
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_VMIX
#undef MY25_MIXING
#define EW_PERIODIC
#define NS_PERIODIC
#define  UV_LDRAG
#define  DIAGNOSTICS_UV
#define  DIAGNOSTICS_TS

#undef BIOLOGY
#ifdef BIOLOGY
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
# define SHORTWAVE
# define ANA_SRFLUX
# define ANA_CLOUD
# define ANA_HUMIDITY
# define ANA_TAIR
#endif

#undef OFFLINE_BIOLOGY
#undef OFFLINE_FLOATS

# elif defined JAVA_1
/*
**  Options for South of Java simulation
*/
 
/* general */

#define CURVGRID
#define MASKING
#define NONLIN_EOS
#define SOLVE3D
#define SALINITY
#define SPLINES
# undef FLOATS
# undef STATIONS
 
/* output stuff */
 
#define NO_WRITE_GRID
#undef OUT_DOUBLE
#define RST_SINGLE
#define AVERAGES
#define AVERAGES_AKT
#define AVERAGES_AKS
#define AVERAGES_AKV
#define AVERAGES_FLUXES
#undef AVERAGES_QUADRATIC
#undef DIAGNOSTICS_TS
#undef DIAGNOSTICS_UV
 
/* advection, dissipation, pressure grad, etc. */
 
#define DJ_GRADPS
 
#define UV_ADV
#define UV_COR
#define UV_SADVECTION
 
#define TS_U3HADVECTION
#define TS_SVADVECTION
 
#define UV_VIS2
#define MIX_S_UV
#define VISC_GRID

          
#define TS_DIF2
#define MIX_GEO_TS
#define DIFF_GRID
 
#define UV_QDRAG
 
/* vertical mixing */
 
#define SOLAR_SOURCE
 
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#  undef LMD_DDMIX
#  undef LMD_BKPP
#  undef LMD_SHAPIRO
#endif
 
# undef GLS_MIXING
# undef MY25_MIXING
 
/* surface forcing */
 
#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE_OUT
# define DIURNAL_SRFLUX
# define EMINUSP
# undef ANA_SRFLUX
# undef ALBEDO
# undef LONGWAVE
#endif
 
/* surface and side corrections */
 
# undef SRELAXATION
# undef QCORRECTION
 
#define SPONGE
# undef TCLIMATOLOGY
# undef TCLM_NUDGING
 
/* Boundary conditions...careful with grid orientation */
 
#undef EASTERN_WALL
#define NORTHERN_WALL
#undef WESTERN_WALL
#undef SOUTHERN_WALL
 
#define RADIATION_2D
 
#ifndef NORTHERN_WALL
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# define NORTH_M3RADIATION
# define NORTH_M3NUDGING
# define NORTH_TRADIATION
# define NORTH_TNUDGING
#endif
 
#ifndef WESTERN_WALL
# define WEST_FSCHAPMAN
# define WEST_M2FLATHER
# define WEST_M3RADIATION
# define WEST_M3NUDGING
# define WEST_TRADIATION
# define WEST_TNUDGING
#endif
 
#ifndef SOUTHERN_WALL
# define SOUTH_FSCHAPMAN
# define SOUTH_M2FLATHER
# define SOUTH_M3RADIATION
# define SOUTH_M3NUDGING
# define SOUTH_TRADIATION
# define SOUTH_TNUDGING
#endif
 
#ifndef EASTERN_WALL
# define EAST_FSCHAPMAN
# define EAST_M2FLATHER
# define EAST_M3RADIATION
# define EAST_M3NUDGING
# define EAST_TRADIATION
# define EAST_TNUDGING
#endif
 
/* roms quirks */
 
#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif
 
# elif defined NEP4
/*
**  Options for Northeast Pacific (NEP4) simulation
*/
 
/* general */

#define CURVGRID
#define MASKING
#define NONLIN_EOS
#define SOLVE3D
#define SALINITY
#define SPLINES
# undef FLOATS
# undef STATIONS
 
/* ice */

#ifdef SOLVE3D
# define  ICE_MODEL
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

/* output stuff */
 
#define NO_WRITE_GRID
#undef OUT_DOUBLE
#define RST_SINGLE
#define AVERAGES
#ifdef SOLVE3D
# define AVERAGES_AKT
# define AVERAGES_AKS
# define AVERAGES_AKV
# define AVERAGES_FLUXES
# undef AVERAGES_QUADRATIC
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
#  undef ALBEDO
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
 
/* point sources (rivers, line sources) */
 
#ifdef SOLVE3D
# define UV_PSOURCE
# define TS_PSOURCE
#endif
 
/* tides */
 
#undef SSH_TIDES
#undef UV_TIDES
#undef ADD_FSOBC
#undef ADD_M2OBC
 
/* Boundary conditions...careful with grid orientation */
 
#if defined NEP4
# define EASTERN_WALL
# define NORTHERN_WALL
#endif
#undef WESTERN_WALL
#undef SOUTHERN_WALL
 
#define RADIATION_2D
 
#ifndef NORTHERN_WALL
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# ifdef SOLVE3D
#  define NORTH_M3RADIATION
#  define NORTH_M3NUDGING
#  define NORTH_TRADIATION
#  define NORTH_TNUDGING
#  define NORTH_MIGRADIENT
# endif
#endif
 
#ifndef WESTERN_WALL
# define WEST_FSCHAPMAN
# define WEST_M2FLATHER
# ifdef SOLVE3D
#  define WEST_M3RADIATION
#  define WEST_M3NUDGING
#  define WEST_TRADIATION
#  define WEST_TNUDGING
#  define WEST_MIGRADIENT
# endif
#endif
 
#ifndef SOUTHERN_WALL
# define SOUTH_FSCHAPMAN
# define SOUTH_M2FLATHER
# ifdef SOLVE3D
#  define SOUTH_M3RADIATION
#  define SOUTH_M3NUDGING
#  define SOUTH_TRADIATION
#  define SOUTH_TNUDGING
#  define SOUTH_MIGRADIENT
# endif
#endif
 
#ifndef EASTERN_WALL
# define EAST_FSCHAPMAN
# define EAST_M2FLATHER
# ifdef SOLVE3D
#  define EAST_M3RADIATION
#  define EAST_M3NUDGING
#  define EAST_TRADIATION
#  define EAST_TNUDGING
#  define EAST_MIGRADIENT
# endif
#endif
 
/* roms quirks */
 
#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif
 
# elif defined NEP_2 || defined NEP_3 || defined NPAC

#if defined NEP_2 || defined NEP_3
# define  ICE_MODEL
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
! NCEP or the other two
# define  NCEP_FLUXES
# undef  ICE_BULK_FLUXES
# undef  CLOUDS
#else
# define SCORRECTION
# define QCORRECTION
# define BULK_FLUXES
# define LONGWAVE
# define ALBEDO
# define ANA_SRFLUX
# define  ICE_BULK_FLUXES
#endif

/*
**  Options for Northeast Pacific Applications.
*/
#define SOLVE3D
#define SPLINES
#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define DJ_GRADPS  /* define if Splines Jacobian pressure gradient */
#undef  TS_A4HADVECTION
#undef  TS_A4VADVECTION
#define TS_U3HADVECTION
#define TS_SVADVECTION
#undef DIAGNOSTICS_UV
#undef DIAGNOSTICS_TS

#define UV_VIS2
#undef  UV_VIS4
#define TS_DIF2
#undef  TS_DIF4
#define VISC_GRID
#define DIFF_GRID
#define MIX_S_UV
#define MIX_GEO_TS

#define NO_WRITE_GRID
#define NONLIN_EOS
#define SALINITY
#define CURVGRID
#define SOLAR_SOURCE
#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS
#define AVERAGES_FLUXES
#define LMD_MIXING
#define LMD_RIMIX
#define LMD_CONVEC
#define LMD_SKPP
#define LMD_NONLOCAL
#undef  MY25_MIXING

#define ANA_BSFLUX
#define ANA_BTFLUX
#undef  ANA_INITIAL
#undef  ANA_SMFLUX
#undef  ANA_STFLUX
#undef  ANA_TCLIMA

#define EASTERN_WALL

# if defined NPAC
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define WESTERN_WALL
#  define MASKING
# else
#  define       SPONGE
#  define       RADIATION_2D
#  undef       TCLIMATOLOGY
#  undef       TCLM_NUDGING
#  undef       M2CLIMATOLOGY
#  undef        M2CLM_NUDGING
#  undef       M3CLIMATOLOGY
#  undef       ZCLIMATOLOGY    /* Processing of SSH climatology */
#  if defined NEP_2 || defined NEP_3
#   define      NORTHERN_WALL
#   define      UV_PSOURCE
#   define      TS_PSOURCE
#   undef       ANA_PSOURCE
#   define      MASKING
#  endif
#  define       SOUTH_FSCHAPMAN
#  define       SOUTH_M2FLATHER
#  define       SOUTH_M3RADIATION
#  define       SOUTH_TRADIATION
#  undef       SOUTH_M2NUDGING
#  define       SOUTH_M3NUDGING
#  define       SOUTH_TNUDGING
#  define       WEST_FSCHAPMAN
#  define       WEST_M2FLATHER
#  define       WEST_M3RADIATION
#  define       WEST_TRADIATION
#  undef       WEST_M2NUDGING
#  define       WEST_M3NUDGING
#  define       WEST_TNUDGING
# endif

# elif defined NP_10

/*
**  Options for North Pacific simulation
*/

#define SOLVE3D
#define SPLINES
#define STATIONS
#define NONLIN_EOS
#define CURVGRID
#define NO_WRITE_GRID
#undef  WRITE_WATER       /* define if only writing water points data */
#define AVERAGES
#undef FLOATS
#define MASKING
#define UV_ADV
#define UV_SVADVECTION
#define UV_COR
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_SVADVECTION

#define VISC_GRID
#define UV_VIS2
#define DIFF_GRID
#define TS_DIF2
#define MIX_S_UV
#define MIX_GEO_TS
#undef SPONGE

#undef UV_PSOURCE
#undef TS_PSOURCE
#undef SSH_TIDES
#undef UV_TIDES
#undef ADD_FSOBC
#undef ADD_M2OBC

#ifdef SOLVE3D

# define SALINITY

# define LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  define LMD_NONLOCAL
#  undef LMD_BKPP
#  undef LMD_SHAPIRO
# endif

# undef GLS_MIXING
# undef MY25_MIXING

# undef TCLIMATOLOGY
# undef TCLM_NUDGING
# define SOLAR_SOURCE
# define DIURNAL_SRFLUX

# define NCEP_FORCE
# ifdef NCEP_FORCE
#  define BULK_FLUXES
#  define LONGWAVE
#  define ALBEDO
#  define ANA_SRFLUX
# else
#  define SRELAXATION
#  define QCORRECTION
# endif

#endif

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX     /* analytical surface momentum stress */
#endif

/* Boundary conditions...careful with grid orientation */

#undef RADIATION_2D

#define EASTERN_WALL
#define NORTHERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL

#ifndef NORTHERN_WALL
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# define NORTH_M3RADIATION
# define NORTH_M3NUDGING
# define NORTH_TRADIATION
# define NORTH_TNUDGING
#endif

#ifndef WESTERN_WALL
# define WEST_FSCHAPMAN
# undef WEST_M2FLATHER
# define WEST_M3RADIATION
# define WEST_M3NUDGING
# define WEST_TRADIATION
# define WEST_TNUDGING
#endif

#ifndef SOUTHERN_WALL
# define SOUTH_FSCHAPMAN
# undef SOUTH_M2FLATHER
# define SOUTH_M3RADIATION
# define SOUTH_M3NUDGING
# define SOUTH_TRADIATION
define SOUTH_TNUDGING
#endif

#ifndef EASTERN_WALL
# define EAST_FSCHAPMAN
# define EAST_M2FLATHER
# define EAST_M3RADIATION
# define EAST_M3NUDGING
# define EAST_TRADIATION
# define EAST_TNUDGING
#endif

#elif defined CGOA
/*
**  Options for Coastal Gulf of Alaska FSL/SMS test problem.
*/

#define UV_ADV
#define UV_VIS2
#define UV_COR
#define UV_PSOURCE

#define TS_U3HADVECTION
#define TS_SVADVECTION
#undef  TS_A4HADVECTION
#undef  TS_A4VADVECTION
#define TS_DIF2
#define SALINITY
#define NONLIN_EOS
#define SOLAR_SOURCE
#define TS_PSOURCE
#define FLOATS

#define DJ_GRADPS
#define SOLVE3D
#define CURVGRID
#define MASKING
#define SPLINES
#define ANA_BSFLUX
#define ANA_BTFLUX

#define MIX_S_UV
#define MIX_GEO_TS
#define LMD_MIXING
#define LMD_RIMIX
#define LMD_CONVEC
#define LMD_SKPP
#undef LMD_BKPP
#define LMD_NONLOCAL

# define NCEP_FORCE
# ifdef NCEP_FORCE
#  define BULK_FLUXES
#  define LONGWAVE
#  define ALBEDO
#  define ANA_SRFLUX
# else
#  define SCORRECTION
#  define QCORRECTION
# endif

#define SPONGE
#define EASTERN_WALL
#define NORTHERN_WALL
#define RADIATION_2D
#define SOUTH_FSCHAPMAN
#define SOUTH_M2FLATHER
#define SOUTH_M3RADIATION
#define SOUTH_TRADIATION
#define SOUTH_M3NUDGING
#define SOUTH_TNUDGING
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_M3RADIATION
#define WEST_TRADIATION
#define WEST_M3NUDGING
#define WEST_TNUDGING

#undef M2CLIMATOLOGY
#undef M3CLIMATOLOGY
#define TCLIMATOLOGY
#undef ZCLIMATOLOGY
#define TCLM_NUDGING
#define NO_WRITE_GRID

# elif defined CGOA_3 || defined CGOA_1 || defined BERING || defined SEBS

/*
**  Options for 3 km Coastal Gulf of Alaska
*/

#define SOLVE3D
#ifdef SOLVE3D
# if defined BERING
#  define  ICE_MODEL
#  undef  ICE_THERMO
#  undef  ICE_MK
#  undef   ICE_ALB_EC92
#  undef   ICE_SMOOTH
#  define  ICE_MOMENTUM
#  define  ICE_MOM_BULK
#  define  ICE_EVP
#  define  ICE_ADVECT
#  define  ICE_SMOLAR
#  define  ICE_UPWIND
!#  define  ICE_BULK_FLUXES
!#  define  BULK_FLUXES
#  define  CLOUDS
# else
#  define SCORRECTION
#  define QCORRECTION
#  define BULK_FLUXES
#  define LONGWAVE
#  define ALBEDO
#  define ANA_SRFLUX
# endif
#else
# define ANA_INITIAL
#endif

#undef OFFLINE_FLOATS

#define SPLINES
#undef FLOATS
#define MASKING
#define UV_ADV
#define UV_COR
#define UV_QDRAG
#ifdef SOLVE3D
# define DJ_GRADPS
# define TS_U3HADVECTION
# define TS_SVADVECTION
# undef  TS_C2HADVECTION
# undef  TS_C2VADVECTION
#endif

!#if defined OFFLINE_FLOATS
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
!# define ANA_VMIX
# define ANA_SSFLUX
!#else
# ifdef SOLVE3D
#  define ADD_FSOBC
#  define ADD_M2OBC
!#  define UV_PSOURCE
!#  define TS_PSOURCE
#  define LMD_MIXING
#  define LMD_RIMIX
#  undef LMD_SHAPIRO
#  undef LMD_CONVEC
#  define LMD_SKPP
#  define LMD_BKPP
#  define LMD_NONLOCAL
#  undef MY25_MIXING
#  define SOLAR_SOURCE
#  undef NCEP_FORCE
#  if defined NCEP_FORCE || defined BULK_FLUXES
#   define LONGWAVE
#   define ALBEDO
#   define ANA_SRFLUX
#  else
#   define SCORRECTION
#   define QCORRECTION
#  endif
# endif
# define SSH_TIDES
# define UV_TIDES
# undef TIDES_ASTRO
!#endif

#define UV_VIS2
#define MIX_S_UV
#define  STATIONS        /* define if writing out station data */
#undef  STATIONS_CGRID  /* define if extracting data at native C-grid */

#define NO_WRITE_GRID
#define NONLIN_EOS
#define CURVGRID
#undef AVERAGES
#ifdef SOLVE3D
# define SPONGE
#endif

#ifdef SOLVE3D
# define SALINITY
# undef TS_DIF2
# undef MIX_GEO_TS
# define FILTERED
# define AVERAGES_AKV
# define AVERAGES_AKT
# define AVERAGES_AKS

# undef TCLIMATOLOGY
# undef TCLM_NUDGING

#endif

#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX     /* analytical surface momentum stress */
# define ANA_INITIAL
#endif

#ifdef OFFLINE_FLOATS
# define EASTERN_WALL
# define NORTHERN_WALL
# define  WESTERN_WALL
# define  SOUTHERN_WALL
#else
# define RADIATION_2D
# if defined CGOA_1 || defined CGOA_3
#  define EASTERN_WALL
#  define NORTHERN_WALL
# endif
# undef  WESTERN_WALL
# undef  SOUTHERN_WALL
#endif

#ifndef EASTERN_WALL
# define EAST_FSCHAPMAN
# define EAST_M2FLATHER
# ifdef SOLVE3D
#  define EAST_M3RADIATION
#  define EAST_M3NUDGING
#  define EAST_TRADIATION
#  define EAST_TNUDGING
#  define EAST_MIGRADIENT
# endif
#endif

#ifndef NORTHERN_WALL
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# ifdef SOLVE3D
#  define NORTH_M3RADIATION
#  define NORTH_M3NUDGING
#  define NORTH_TRADIATION
#  define NORTH_TNUDGING
#  define NORTH_MIGRADIENT
# endif
#endif

#ifndef SOUTHERN_WALL
# define SOUTH_FSCHAPMAN
# define SOUTH_M2FLATHER
# ifdef SOLVE3D
#  define SOUTH_M3RADIATION
#  define SOUTH_M3NUDGING
#  define SOUTH_TRADIATION
#  define SOUTH_TNUDGING
#  define SOUTH_MIGRADIENT
# endif
#endif

#ifndef WESTERN_WALL
# define WEST_FSCHAPMAN
# define WEST_M2FLATHER
# ifdef SOLVE3D
#  define WEST_M3RADIATION
#  define WEST_M3NUDGING
#  define WEST_TRADIATION
#  define WEST_TNUDGING
#  define WEST_MIGRADIENT
# endif
#endif

# elif defined WBC_1 || defined WBC_2 || defined WBC_3

/*
**  Options for Western Boundary Current example:
*/

#define UV_VIS2
#define UV_COR
#undef AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX

#  if defined WBC_2 || defined WBC_3
#define MASKING
#define ANA_MASK
#  endif

# elif defined CRIT_LAT

/*
**  Options for critical latitude test.
*/

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define SPONGE
#define SOLVE3D
#define SPLINES
#undef AVERAGES

#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_M3OBC

#define EW_PERIODIC
#undef  RADIATION_2D
#define NORTH_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_M3RADIATION
#define SOUTH_FSCLAMPED
#define SOUTH_M2CLAMPED
#define SOUTH_M3CLAMPED

#endif
