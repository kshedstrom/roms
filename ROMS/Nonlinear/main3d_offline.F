#include "cppdefs.h"
#if defined NONLINEAR && defined SOLVE3D && \
    (defined OFFLINE_BIOLOGY || defined OFFLINE_FLOATS)
      SUBROUTINE main3d_offline (ng)
!
!svn $Id: main3d_offline.F 982 2009-05-22 20:04:01Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine is the main driver for nonlinear ROMS/TOMS when     !
!  configurated as a full 3D baroclinic ocean model.  It  advances     !
!  forward the primitive equations for a single time step.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef MODEL_COUPLING
      USE mod_coupler
# endif
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
# ifdef ANA_VMIX
      USE analytical_mod, ONLY : ana_vmix
# endif
# ifdef BIOLOGY
      USE biology_mod, ONLY : biology
# endif
# ifdef BBL_MODEL
      USE bbl_mod, ONLY : bblm
# endif
# ifdef BULK_FLUXES
#  ifdef CCSM_FLUXES
      USE ccsm_flux_mod, ONLY : ccsm_flux
#  else
      USE bulk_flux_mod, ONLY : bulk_flux
#  endif
# endif
# if defined ALBEDO || defined NCEP_FLUXES
      USE cawdir_eval_mod, ONLY : cawdir_eval
# endif
# if defined NCEP_FLUXES
      USE ncep_flux_mod, ONLY : ncep_flux
# endif
# ifdef BVF_MIXING
      USE bvf_mix_mod, ONLY : bvf_mix
# endif
      USE diag_mod, ONLY : diag
# ifdef TLM_CHECK
      USE dotproduct_mod, ONLY : nl_dotproduct
# endif
# if defined W4DPSAS || defined NLM_OUTER
      USE forcing_mod, ONLY : forcing
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE frc_adjust_mod, ONLY : frc_NLadjust, load_NLforcing
# endif
# ifdef GLS_MIXING
      USE gls_corstep_mod, ONLY : gls_corstep
      USE gls_prestep_mod, ONLY : gls_prestep
# endif
# if defined DIFF_3DCOEF || defined VISC_3DCOEF
      USE hmixing_mod, ONLY : hmixing
# endif
      USE ini_fields_mod, ONLY : ini_fields, ini_zeta
# ifdef LMD_MIXING
      USE lmd_vmix_mod, ONLY : lmd_vmix
# endif
# ifdef MY25_MIXING
      USE my25_corstep_mod, ONLY : my25_corstep
      USE my25_prestep_mod, ONLY : my25_prestep
# endif
# ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2atm_coupling
# endif
# ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2wav_coupling
# endif
# ifdef ASSIMILATION
      USE oi_update_mod, ONLY : oi_update
# endif
      USE omega_mod, ONLY : omega
# ifdef NEARSHORE_MELLOR
      USE radiation_stress_mod, ONLY : radiation_stress
# endif
# ifndef TS_FIXED
      USE rho_eos_mod, ONLY : rho_eos
# endif
      USE rhs3d_mod, ONLY : rhs3d
# ifdef SEDIMENT
      USE sediment_mod, ONLY : sediment
# endif
# if defined AVERAGES && !defined ADJOINT
      USE set_avg_mod, ONLY : set_avg
# endif
# if defined AVERAGES2 && !defined ADJOINT
      USE set_avg2_mod, ONLY : set_avg2
# endif
# if defined ICE_MODEL && defined ICE_THERMO
      USE ice_frazil_mod, ONLY : ice_frazil
# endif
      USE set_massflux_mod, ONLY : set_massflux
# if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
      USE set_tides_mod, ONLY : set_tides
# endif
      USE set_vbc_mod, ONLY : set_vbc
# ifdef SET_ZETA
      USE set_zeta_mod, ONLY : set_zeta
# endif
      USE step2d_mod, ONLY : step2d
# ifndef TS_FIXED
      USE step3d_t_mod, ONLY : step3d_t
# endif
      USE step3d_uv_mod, ONLY : step3d_uv
# ifdef FLOATS
      USE step_floats_mod, ONLY : step_floats
# endif
      USE wvelocity_mod, ONLY : wvelocity
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: my_iif, next_indx1, subs, tile, thread
# if defined FLOATS
      integer :: Lend, Lstr, chunk_size
# endif
!
!=======================================================================
!  Time-step nonlinear 3D primitive equations.
!=======================================================================
!
!  Set time indices and time clock.
!
      nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
      nnew(ng)=3-nstp(ng)
      nrhs(ng)=nstp(ng)
      CALL time_string (time(ng), time_code(ng))
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL set_data (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
      IF (exit_flag.ne.NoError) RETURN

# if defined W4DPSAS || defined NLM_OUTER
!
!-----------------------------------------------------------------------
!  If appropriate, add convolved adjoint solution impulse forcing to
!  the nonlinear model solution. Notice that the forcing is only needed
!  after finishing all the inner loops. The forcing is continuous.
!  That is, it is time interpolated at every time-step from available
!  snapshots (FrequentImpulse=TRUE).
!-----------------------------------------------------------------------
!
      IF (FrequentImpulse) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL forcing (ng, TILE, kstp(ng), nstp(ng))
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Initialize all time levels and compute other initial fields.
!-----------------------------------------------------------------------
!
      IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL ini_zeta (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Initialize other state variables.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL ini_fields (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
!  quatities and report global diagnostics.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL set_massflux (ng, TILE)
# ifndef TS_FIXED
          CALL rho_eos (ng, TILE)
# endif
          CALL diag (ng, TILE)
# ifdef TLM_CHECK
          CALL nl_dotproduct (ng, TILE, Lnew(ng))
# endif
        END DO
      END DO
!$OMP END PARALLEL DO
      IF (exit_flag.ne.NoError) RETURN

# ifdef AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Couple ocean to atmosphere model every "CoupleSteps(Iatmos)"
!  timesteps: get air/sea fluxes.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).ne.ntstart(ng)).and.                                 &
     &    MOD(iic(ng)-1,CoupleSteps(Iatmos,ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL ocn2atm_coupling (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# ifdef WAVES_OCEAN
!
!-----------------------------------------------------------------------
!  Couple to ocean to waves model every "CoupleSteps(Iwaves)"
!  timesteps: get waves/ocean fluxes.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).ne.ntstart(ng)).and.                                 &
     &    MOD(iic(ng)-1,CoupleSteps(Iwaves,ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL ocn2wav_coupling (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# ifdef NEARSHORE_MELLOR
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL radiation_stress (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif
!
# ifndef OFFLINE_FLOATS
!-----------------------------------------------------------------------
!  Set fields for vertical boundary conditions. Process tidal forcing,
!  if any.
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!  Interpolate surface forcing increments and adjust surface forcing.
!  Load surface forcing into storage arrays.
# endif
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
# if defined NCEP_FLUXES || defined ALBEDO
          CALL cawdir_eval(ng, TILE)
# endif
# ifdef BULK_FLUXES
#  ifdef CCSM_FLUXES
          CALL ccsm_flux (ng, TILE)
#  else
          CALL bulk_flux (ng, TILE)
#  endif
# endif
# ifdef NCEP_FLUXES
          CALL ncep_flux(ng, TILE)
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
          CALL frc_NLadjust (ng, TILE, Lfinp(ng), Lfout(ng))
          CALL load_NLforcing (ng, TILE, Lfout(ng))
# endif
# ifdef BBL_MODEL
          CALL bblm (ng, TILE)
# endif
          CALL set_vbc (ng, TILE)
# if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
          CALL set_tides (ng, TILE)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Run ice model for one step
!-----------------------------------------------------------------------
!
#  if defined ICE_MODEL
          CALL seaice(ng)
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute time-dependent vertical/horizontal mixing coefficients for
!  momentum and tracers. Compute S-coordinate vertical velocity,
!  diagnostically from horizontal mass divergence.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,nstp,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
# ifndef OFFLINE_FLOATS
#  if defined ANA_VMIX
          CALL ana_vmix (ng, TILE, iNLM)
#  elif defined LMD_MIXING
          CALL lmd_vmix (ng, TILE)
#  elif defined BVF_MIXING
          CALL bvf_mix (ng, TILE)
#  endif
# endif
# if defined DIFF_3DCOEF || defined VISC_3DCOEF
          CALL hmixing (ng, TILE)
# endif
# if !defined OCLIMATOLOGY
          CALL omega (ng, TILE)
# endif
          CALL wvelocity (ng, TILE, nstp(ng))
!
!--------------------------------------------------------------------
!  Compute deformation-dependent horizontal eddy viscosities
!  using Smagorinsky(1969) formulation
!--------------------------------------------------------------------
!
        END DO
      END DO
!$OMP END PARALLEL DO

# ifdef SET_ZETA
!
!-----------------------------------------------------------------------
!  Set free-surface to it time-averaged value.  If applicable,
!  accumulate time-averaged output data which needs a irreversible
!  loop in shared-memory jobs.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1     ! irreversible loop
          CALL set_zeta (ng, TILE)
#  ifdef DIAGNOSTICS
          CALL set_diags (ng, TILE)
#  endif
#  if defined AVERAGES && !defined ADJOINT
          CALL set_avg (ng, TILE)
#  endif
#  if defined AVERAGES2 && !defined ADJOINT
          CALL set_avg2 (ng, TILE)
#  endif
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
      CALL output (ng)
      IF ((exit_flag.ne.NoError).or.(iic(ng).eq.(ntend(ng)+1))) RETURN
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL rhs3d (ng, TILE)
#  ifdef MY25_MIXING
          CALL my25_prestep (ng, TILE)
#  elif defined GLS_MIXING
          CALL gls_prestep (ng, TILE)
#  endif
        END DO
      END DO
!$OMP END PARALLEL DO

# else
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.  If applicable,
!  accumulate time-averaged output data.
!-----------------------------------------------------------------------
!
# ifndef OFFLINE_FLOATS
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
#  ifdef DIAGNOSTICS
          CALL set_diags (ng, TILE)
#  endif
          CALL rhs3d (ng, TILE)
#  ifdef MY25_MIXING
          CALL my25_prestep (ng, TILE)
#  elif defined GLS_MIXING
          CALL gls_prestep (ng, TILE)
#  endif
        END DO
      END DO
!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
      CALL output (ng)
      IF ((exit_flag.ne.NoError).or.(iic(ng).eq.(ntend(ng)+1))) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Time-step vertical mixing turbulent equations and passive tracer
!  source and sink terms, if applicable.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,nnew,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL omega (ng, TILE)
# ifdef MY25_MIXING
          CALL my25_corstep (ng, TILE)
# elif defined GLS_MIXING
          CALL gls_corstep (ng, TILE)
# endif
# ifdef BIOLOGY
          CALL biology (ng, TILE)
# endif
# ifdef SEDIMENT
          CALL sediment (ng, TILE)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO

# ifndef TS_FIXED
!
!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!
# if !defined OFFLINE_FLOATS
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL step3d_t (ng, TILE)
#  if defined ICE_MODEL && defined ICE_THERMO
          CALL ice_frazil(ng, TILE)
#  endif
        END DO
      END DO
!$OMP END PARALLEL DO
# endif
# endif

# ifdef ASSIMILATION
!
!-----------------------------------------------------------------------
!  Assimilate observations via Optimal Interpolation.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL oi_update (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif

# ifdef FLOATS
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,chunk_size,Lstr,Lend)                  &
!$OMP&            SHARED(ng,numthreads,Nfloats)
      DO thread=0,numthreads-1
        chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
        Lstr=1+thread*chunk_size
        Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
        CALL step_floats (ng, Lstr, Lend)
      END DO
!$OMP END PARALLEL DO
# endif
# if defined FLOATS
!
!  Shift floats time indices.
!
      nfp1(ng)=MOD(nfp1(ng)+1,NFT+1)
      nf(ng)  =MOD(nf(ng)  +1,NFT+1)
      nfm1(ng)=MOD(nfm1(ng)+1,NFT+1)
      nfm2(ng)=MOD(nfm2(ng)+1,NFT+1)
      nfm3(ng)=MOD(nfm3(ng)+1,NFT+1)
# endif
      RETURN
      END SUBROUTINE main3d_offline
#else
      SUBROUTINE main3d_offline
      RETURN
      END SUBROUTINE main3d_offline
#endif
