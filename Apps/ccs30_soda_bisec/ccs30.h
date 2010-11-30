#define OUT_DOUBLE
#define CURVGRID
#define MASKING
#define SOLVE3D
#define SPLINES

#define AVERAGES
#define AVERAGES_FLUXES

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define UV_QDRAG
#define DJ_GRADPS

#define SALINITY
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define NONLIN_EOS

#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE_OUT
# define EMINUSP
# define SPEC_HUM
#else
# define QCORRECTION
# undef SCORRECTION
#endif

#define LMD_MIXING
#ifdef LMD_MIXING
# undef  DIURNAL_SRFLUX
# define SOLAR_SOURCE
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL
# define LMD_SKPP
# undef  LMD_BKPP
#endif

#undef  CLOSED_OBC
#ifdef CLOSED_OBC
# define NORTHERN_WALL
# define SOUTHERN_WALL
# define EASTERN_WALL
# define WESTERN_WALL
#else
# define EASTERN_WALL
# define RADIATION_2D
# define WEST_FSCHAPMAN
# define WEST_M2FLATHER
# define WEST_M3CLAMPED
# define WEST_TCLAMPED
# undef  WEST_M3RADIATION
# undef  WEST_TRADIATION
# undef  WEST_M3NUDGING
# undef  WEST_TNUDGING
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# define NORTH_M3CLAMPED
# define NORTH_TCLAMPED
# undef  NORTH_M3RADIATION
# undef  NORTH_TRADIATION
# undef  NORTH_M3NUDGING
# undef  NORTH_TNUDGING
# define SOUTH_FSCHAPMAN
# define SOUTH_M2FLATHER
# define SOUTH_M3CLAMPED
# define SOUTH_TCLAMPED
# undef  SOUTH_M3RADIATION
# undef  SOUTH_TRADIATION
# undef  SOUTH_M3NUDGING
# undef  SOUTH_TNUDGING
#endif

#define SPONGE
#ifdef SPONGE
# define CCS_SPONGE
#endif

#define NEMURO
#ifdef NEMURO
# define NEMURO_SAN
# ifdef NEMURO_SAN
#  undef  FISH_FEEDBACK
#  define PREDATOR
#  define FLEET
#  undef  ANA_SPAWN_DIST
#  define EGGS_BISECTION
#  undef  EGGS_TREE_FORT
#  undef  EGGS_TREE_CXX
#  undef  EGGS_VECTOR_CXX
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
