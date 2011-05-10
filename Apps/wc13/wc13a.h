#define OUT_DOUBLE

#undef  FORWARD_WRITE
#undef  AD_SENSITIVITY

#ifdef  AD_SENSITIVITY
# define FORWARD_READ
# define FORWARD_MIXING
#endif

#define AVERAGES
#define AVERAGES_FLUXES

#undef DIAGNOSTICS_UV
#undef DIAGNOSTICS_TS
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define UV_U3HADVECTION
#define UV_SADVECTION
#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define TS_DIF2
#define MIX_GEO_TS
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define MASKING
#define SOLVE3D
#define SPLINES

#define UV_QDRAG
#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE_OUT
# define EMINUSP
# define SPEC_HUM
#else
# define QCORRECTION
# undef SCORRECTION
#endif

#undef  MY25_MIXING
#define LMD_MIXING
#ifdef  AD_SENSITIVITY
# undef  LMD_MIXING
#endif
#ifdef LMD_MIXING
# undef  DIURNAL_SRFLUX
# define SOLAR_SOURCE
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL
# define LMD_SKPP
#endif

#undef  M2CLIMATOLOGY
#undef  M3CLIMATOLOGY
#undef  TCLIMATOLOGY
#undef  M2CLM_NUDGING
#undef  M3CLM_NUDGING
#undef  TCLM_NUDGING

#undef CLOSED_OBC
#ifdef CLOSED_OBC
# define NORTHERN_WALL
# define SOUTHERN_WALL
# define EASTERN_WALL
# define WESTERN_WALL
#else
# define EASTERN_WALL
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# define NORTH_M3CLAMPED
# define NORTH_TCLAMPED
# define SOUTH_FSCHAPMAN
# define SOUTH_M2FLATHER
# define SOUTH_M3CLAMPED
# define SOUTH_TCLAMPED
# define WEST_FSCHAPMAN
# define WEST_M2FLATHER
# define WEST_M3CLAMPED
# define WEST_TCLAMPED
#endif

#undef  FLOATS

#define NEMURO
#ifdef NEMURO
# define NEMURO_SAN
# ifdef NEMURO_SAN
#  define FISH_FEEDBACK
#  define PREDATOR
#  undef FISHING_FLEET
#  define ANA_SPAWN_DIST
#  define EGGS_BISECTION
#  undef EGGS_TREE_FORT
#  undef EGGS_TREE_CXX
#  undef EGGS_VECTOR_CXX
# endif
# define BIO_SEDIMENT
# define HOLLING_GRAZING
# undef  IVLEV_EXPLICIT
# undef  IRON_LIMIT
# undef  IRON_RELAX
# undef  IRON_RSIN
# undef  ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
# undef  ANA_SRFLUX
# undef  DIAGNOSTICS_BIO
#endif

#define ANA_BSFLUX
#define ANA_BTFLUX
