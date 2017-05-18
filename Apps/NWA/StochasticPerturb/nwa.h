/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2013 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for NWA simulation
*/

#undef NETCDF4
#undef PARALLEL_IO
#undef OFFLINE_FLOATS

/* general */

#define CURVGRID
#define MASKING
#define NONLIN_EOS
#define SOLVE3D
#define SALINITY
#ifdef SOLVE3D
/* #define SPLINES */
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#endif
#undef FLOATS
#undef STATIONS
#undef WET_DRY

#undef T_PASSIVE
#ifdef T_PASSIVE
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
# define ANA_PASSIVE
# define TRC_PSOURCE
# define ANA_TRC_PSOURCE
# define AGE_PASSIVE
#endif

/* output stuff */

 
#define NO_HIS
#define LONG_NUMS
#define NO_WRITE_GRID
#undef OUT_DOUBLE
#define RST_SINGLE
#define AVERAGES
#undef AVERAGES2
#ifdef SOLVE3D
# undef AVERAGES_DETIDE
# undef DIAGNOSTICS_TS
#endif
#undef DIAGNOSTICS_UV

/* advection, dissipation, pressure grad, etc. */

#ifdef SOLVE3D
# define DJ_GRADPS               /* use splines density Jacobian (Shchepetkin, 2000) for pressure gradient algorithm */
#endif

#define UV_ADV
#define UV_COR
#undef UV_SADVECTION             /* deep water so no splines vertical advection */
#define UV_VIS2                  /* horizontal, harmonic viscosity of momentum */
#undef UV_SMAGORINSKY           /* use Smagorinsky-like viscosity */
#define VISC_3DCOEF
#define MIX_S_UV                 /* mixing along constant S-surfaces for horizontal mixing of momentum */
#define VISC_GRID
#undef SPONGE                   /* impose a sponge layer near the lateral boundary */

#ifdef SOLVE3D
# define TS_DIF2
# define MIX_GEO_TS              /* mixing on geopotential surfaces for horizontal mixing of tracers */
# define DIFF_GRID               /* scale diffusion coefficients by grid size for horizontal mixing of tracers */
# define TS_U3HADVECTION         /* 3rd-order upstream biased advection */
# define TS_C4VADVECTION         /* 4th-order center vert advection instead */
# undef TS_MPDATA
#endif

/* vertical mixing */

#ifdef SOLVE3D
# define LMD_MIXING              /* activate Large et al. (1994) interior closure for vertical turbulent mixing of momentum and tracers */
# ifdef LMD_MIXING
#  define LMD_RIMIX              /* add diffusivity due to shear instability for the Large et al. (1994) K-profile parameterization mixing */
#  define LMD_CONVEC             /* add convective mixing due to shear instability for the Large et al. (1994) K-profile parameterization mixing */
#  define LMD_SKPP               /* use surface boundary layer KPP mixing for the Large et al. (1994) K-profile parameterization mixing */
#  define LMD_BKPP
#  define LMD_NONLOCAL           /* use nonlocal transport for the Large et al. (1994) K-profile parameterization mixing */
#  define LMD_SHAPIRO            /* use Shapiro filtering boundary layer depth for the Large et al. (1994) K-profile parameterization mixing */
#  undef LMD_DDMIX
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
# define CORE_FORCING
# define BULK_FLUXES
# define CCSM_FLUXES
# if defined BULK_FLUXES || defined CCSM_FLUXES
#  define LONGWAVE_OUT
#  undef DIURNAL_SRFLUX
#  define EMINUSP
#  undef ANA_SRFLUX
#  undef ALBEDO
#  define SOLAR_SOURCE
#  define ALBEDO_CURVE
#  undef ALBEDO_FILE
#  undef LONGWAVE
# endif
#endif
# define WTYPE_GRID

/* surface and side corrections */

#ifdef SOLVE3D
# define SCORRECTION /* climate run : no SSS restoring */
# define SSSC_THRESHOLD
# undef SRELAXATION
# undef QCORRECTION
#endif

/* Limitation of temperature */
#define TEMP_CLAMP
 
/* commented for conservation test */
#ifdef SOLVE3D
# undef ANA_TCLIMA
# undef ANA_NUDGCOEF
#endif

/* tides */


#define LTIDES
#ifdef LTIDES
# define FILTERED
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# undef RAMP_TIDES
# define TIDES_ASTRO
# define POT_TIDES
# define UV_LDRAG
# define UV_DRAG_GRID
# define ANA_DRAG
# define LIMIT_BSTRESS
# undef UV_QDRAG
#else
# define UV_QDRAG
# undef M2TIDE_DIFF
#endif

/* Stochastic perturbations */
#define STOPERTURB
#ifdef STOPERTURB
# define FORCE_PERTURB
# define EOS_PERTURB
#endif


/* point sources (rivers, line sources) */

/* Using Runoff instead now */
#ifdef SOLVE3D
# define RUNOFF
# undef UV_PSOURCE
# undef ANA_PSOURCE
# undef TS_PSOURCE
#endif

#define RADIATION_2D


/* roms quirks */

#ifdef SOLVE3D
# define ANA_BSFLUX             /* set kinematic bottom salinity flux with an analytical expression */
# define ANA_BTFLUX             /* set kinematic bottom heat flux with an analytical expression */
#else
# define ANA_SMFLUX
#endif

/* biological model options */

#undef BIO_UMAINE       /* Chai et al. (2002) CoSINE model */

#if defined BIO_UMAINE
# define OPTIC_UMAINE
# define OXYGEN
# define CARBON
# undef ANA_BIOLOGY        /* analytical biology initial conditions */
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
#endif
