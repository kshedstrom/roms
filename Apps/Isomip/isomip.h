/*
**  Ice Shelf test 2.01 as defined by J R Hunter
*/

#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#define UV_VIS2
#define UV_QDRAG
#undef  MIX_GEO_UV
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#undef  MIX_S_TS
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define SPHERICAL
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#define ICESHELF
#undef  AVERAGES

#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
#undef  ANA_VMIX
#define LMD_MIXING
#ifdef  LMD_MIXING
# undef  LMD_RIMIX
# define LMD_CONVEC
# undef  LMD_SKPP
# undef  LMD_NONLOCAL
#endif
