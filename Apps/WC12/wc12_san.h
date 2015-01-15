#define WC12

#define OUT_DOUBLE
#define AVERAGES

#define CURVGRID
#define MASKING
#define SOLVE3D
#define SPLINES
#define DJ_GRADPS

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define UV_QDRAG
#undef DIAGNOSTICS_UV
#undef  TS_MPDATA
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define SALINITY
#define NONLIN_EOS
#undef DIAGNOSTICS_TS
#define ANA_BTFLUX
#define ANA_BSFLUX

#define SPONGE
#define ANA_NUDGCOEF

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
# define LMD_BKPP
# undef  LMD_SHAPIRO
# undef  RI_HORAVG
# undef  RI_VERAVG
#endif



#define NEMURO
#ifdef NEMURO
# define NEMURO_SAN
# ifdef NEMURO_SAN
#  define FISH_FEEDBACK
#  define PREDATOR
#  define FISHING_FLEET
#  undef  ANA_SPAWN_DIST
#  define EGGS_BISECTION
#  undef  EGGS_TREE_FORT
#  undef  EGGS_TREE_CXX
#  undef  EGGS_VECTOR_CXX
# endif
# undef  BIO_SEDIMENT
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

