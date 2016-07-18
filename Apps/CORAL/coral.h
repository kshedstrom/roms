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
**  Options for CORAL simulation
*/

#define NO_HIS
#undef NETCDF4
#undef PARALLEL_IO
#undef OFFLINE_FLOATS

/* general */

#define CURVGRID
#define MASKING
#define NONLIN_EOS
#define SOLVE3D
#define SALINITY
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES
#undef FLOATS
#define STATIONS
#undef WET_DRY

/* output stuff */

#define NO_WRITE_GRID
#undef OUT_DOUBLE
#define RST_SINGLE
#define AVERAGES
#undef AVERAGES2
#undef AVERAGES_DETIDE
#undef DIAGNOSTICS_TS
#undef DIAGNOSTICS_UV

/* advection, dissipation, pressure grad, etc. */
#ifdef SOLVE3D
# define DJ_GRADPS               /* use splines density Jacobian (Shchepetkin, 2000) for pressure gradient algorithm */
#endif

#define UV_ADV
#define UV_COR
#undef UV_SADVECTION             /* deep water so no splines vertical advection */
#define UV_VIS2                  /* horizontal, harmonic viscosity of momentum */
#undef UV_SMAGORINSKY            /* use Smagorinsky-like viscosity */
#define VISC_3DCOEF
#define MIX_S_UV                 /* mixing along constant S-surfaces for horizontal mixing of momentum */
#define VISC_GRID
#define SPONGE                   /* impose a sponge layer near the lateral boundary */

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
# define SOLAR_SOURCE
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
#  define DIURNAL_SRFLUX
#  define EMINUSP
#  undef ANA_SRFLUX
#  undef ALBEDO_CLOUD
#  define ALBEDO_CURVE
#  undef LONGWAVE
# endif
#endif

/* surface and side corrections */

#ifdef SOLVE3D
# define SCORRECTION
# define SSSC_THRESHOLD
# undef QCORRECTION
#endif

#ifdef SOLVE3D
# undef TCLIMATOLOGY
# undef TCLM_NUDGING
#endif

/* point sources (rivers, line sources) */
/* Using Runoff instead now */

#ifdef SOLVE3D
# define RUNOFF
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
#endif

#define RADIATION_2D

/* roms quirks */

#ifdef SOLVE3D
# define ANA_BSFLUX             /* set kinematic bottom salinity flux with an analytical expression */
# define ANA_BTFLUX             /* set kinematic bottom heat flux with an analytical expression */
#else
# define ANA_SMFLUX
#endif
