#include "cppdefs.h"
      MODULE set_avg2_mod
#if defined AVERAGES2 && (!defined ADJOINT && defined NONLINEAR)
!
!svn $Id: set_avg.F 702 2008-08-12 16:44:47Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beggining of the next time-step.                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: set_avg2

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_avg2 (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_average2
# if defined FORWARD_WRITE && defined SOLVE3D
      USE mod_coupling
# endif
      USE mod_forces
# ifdef SOLVE3D
      USE mod_grid
# endif
# ifdef ICE_MODEL
      USE mod_ice
# endif
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 5)
# endif
      CALL set_avg2_tile (ng, tile,                                     &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   KOUT,                                          &
# ifdef SOLVE3D
#  ifdef ICE_MODEL
     &                   IOUT, IUOUT, IEOUT,                            &
#  endif
     &                   NOUT,                                          &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
     &                   OCEAN(ng) % t,                                 &
#  ifdef BIO_GOANPZ
     &                   OCEAN(ng) % st,                                &
#  endif
     &                   OCEAN(ng) % rho,                               &
#  ifdef LMD_SKPP
     &                   MIXING(ng) % hsbl,                             &
#  endif
#  ifdef LMD_BKPP
     &                   MIXING(ng) % hbbl,                             &
#  endif
#  ifdef ICE_MODEL
     &                   ICE(ng)%ui,                                    &
     &                   ICE(ng)%vi,                                    &
     &                   ICE(ng)%ai,                                    &
     &                   ICE(ng)%hi,                                    &
     &                   ICE(ng)%hsn,                                   &
     &                   ICE(ng)%tis,                                   &
     &                   ICE(ng)%ti,                                    &
     &                   ICE(ng)%sfwat,                                 &
     &                   ICE(ng)%ageice,                                &
     &                   ICE(ng)%io_mflux,                              &
     &                   ICE(ng)%sig11,                                 &
     &                   ICE(ng)%sig12,                                 &
     &                   ICE(ng)%sig22,                                 &
     &                   ICE(ng)%t0mk,                                  &
     &                   ICE(ng)%s0mk,                                  &
     &                   ICE(ng)%utau_iw,                               &
     &                   ICE(ng)%chu_iw,                                &
#  endif
#  ifdef AVERAGES_FLUXES
     &                   FORCES(ng) % stflx,                            &
#   ifdef BULK_FLUXES
     &                   FORCES(ng) % lhflx,                            &
     &                   FORCES(ng) % shflx,                            &
     &                   FORCES(ng) % lrflx,                            &
     &                   FORCES(ng) % Uwind,                            &
     &                   FORCES(ng) % Vwind,                            &
#    ifdef EMINUSP
     &                   FORCES(ng) % evap,                             &
     &                   FORCES(ng) % rain,                             &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                   FORCES(ng) % srflx,                            &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr,                            &
# endif
     &                   OCEAN(ng) % ubar,                              &
     &                   OCEAN(ng) % vbar,                              &
     &                   OCEAN(ng) % zeta,                              &
# ifdef SOLVE3D
     &                   AVERAGE2(ng) % avgu3d,                          &
     &                   AVERAGE2(ng) % avgv3d,                          &
     &                   AVERAGE2(ng) % avgt,                            &
#  ifdef BIO_GOANPZ
     &                   AVERAGE2(ng) % avgst,                           &
#  endif
     &                   AVERAGE2(ng) % avgrho,                          &
#  ifdef LMD_SKPP
     &                   AVERAGE2(ng) % avghsbl,                         &
#  endif
#  ifdef LMD_BKPP
     &                   AVERAGE2(ng) % avghbbl,                         &
#  endif
#  ifdef ICE_MODEL
     &                   AVERAGE2(ng) % avguice,                         &
     &                   AVERAGE2(ng) % avgvice,                         &
     &                   AVERAGE2(ng) % avgaice,                         &
     &                   AVERAGE2(ng) % avghice,                         &
     &                   AVERAGE2(ng) % avghsno,                         &
     &                   AVERAGE2(ng) % avgtice,                         &
     &                   AVERAGE2(ng) % avgtimid,                        &
     &                   AVERAGE2(ng) % avgsfwat,                        &
     &                   AVERAGE2(ng) % avgageice,                       &
     &                   AVERAGE2(ng) % avgiomflx,                       &
     &                   AVERAGE2(ng) % avgsig11,                        &
     &                   AVERAGE2(ng) % avgsig12,                        &
     &                   AVERAGE2(ng) % avgsig22,                        &
     &                   AVERAGE2(ng) % avgT0mk,                         &
     &                   AVERAGE2(ng) % avgS0mk,                         &
     &                   AVERAGE2(ng) % avgutau_iw,                      &
     &                   AVERAGE2(ng) % avgchu_iw,                       &
#  endif
#  ifdef AVERAGES_FLUXES
     &                   AVERAGE2(ng) % avgstf,                          &
     &                   AVERAGE2(ng) % avgswf,                          &
#   ifdef BULK_FLUXES
     &                   AVERAGE2(ng) % avglhf,                          &
     &                   AVERAGE2(ng) % avgshf,                          &
     &                   AVERAGE2(ng) % avglrf,                          &
     &                   AVERAGE2(ng) % avguwind,                        &
     &                   AVERAGE2(ng) % avgvwind,                        &
#    ifdef EMINUSP
     &                   AVERAGE2(ng) % avgevap,                         &
     &                   AVERAGE2(ng) % avgrain,                         &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                   AVERAGE2(ng) % avgsrf,                          &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                   AVERAGE2(ng) % avgsus,                          &
     &                   AVERAGE2(ng) % avgsvs,                          &
     &                   AVERAGE2(ng) % avgbus,                          &
     &                   AVERAGE2(ng) % avgbvs,                          &
# endif
     &                   AVERAGE2(ng) % avgu2d,                          &
     &                   AVERAGE2(ng) % avgv2d,                          &
     &                   AVERAGE2(ng) % avgzeta)

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_avg2
!
!***********************************************************************
      SUBROUTINE set_avg2_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Kout,                                    &
# ifdef SOLVE3D
#  ifdef ICE_MODEL
     &                         Iout, Iuout, Ieout,                      &
#  endif
     &                         Nout,                                    &
     &                         pm, pn,                                  &
     &                         u, v,                                    &
     &                         t,                                       &
# ifdef BIO_GOANPZ
     &                         st,                                      &
# endif
     &                         rho,                                     &
#  ifdef LMD_SKPP
     &                         hsbl,                                    &
#  endif
#  ifdef LMD_BKPP
     &                         hbbl,                                    &
#  endif
#  ifdef ICE_MODEL
     &                         ui, vi, ai, hi, hsn, tis, ti, sfwat,     &
     &                         ageice, io_mflux, sig11, sig12, sig22,   &
     &                         t0mk, s0mk, utau_iw, chu_iw,             &
#  endif
#  ifdef AVERAGES_FLUXES
     &                         stflx,                                   &
#   ifdef BULK_FLUXES
     &                         lhflx, shflx, lrflx, Uwind, Vwind,       &
#    ifdef EMINUSP
     &                         evap, rain,                              &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                         srflx,                                   &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                         sustr, svstr, bustr, bvstr,              &
# endif
     &                         ubar, vbar,                              &
     &                         zeta,                                    &
# ifdef SOLVE3D
     &                         avgu3d, avgv3d,                          &
     &                         avgt,                                    &
# ifdef BIO_GOANPZ
     &                         avgst,                                   &
# endif
     &                         avgrho,                                  &
#  ifdef LMD_SKPP
     &                         avghsbl,                                 &
#  endif
#  ifdef LMD_BKPP
     &                         avghbbl,                                 &
#  endif
#  ifdef ICE_MODEL
     &                         avguice, avgvice, avgaice, avghice,      &
     &                         avghsno, avgtice, avgtimid, avgsfwat,    &
     &                         avgageice, avgiomflx, avgsig11,          &
     &                         avgsig12, avgsig22, avgT0mk, avgS0mk,    &
     &                         avgutau_iw, avgchu_iw,                   &
#  endif
#  ifdef AVERAGES_FLUXES
     &                         avgstf, avgswf,                          &
#   ifdef BULK_FLUXES
     &                         avglhf, avgshf, avglrf,                  &
     &                         avguwind, avgvwind,                      &
#    ifdef EMINUSP
     &                         avgevap, avgrain,                        &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                         avgsrf,                                  &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                         avgsus, avgsvs, avgbus, avgbvs,          &
# endif
     &                         avgu2d, avgv2d,                          &
     &                         avgzeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kout
# ifdef SOLVE3D
#  ifdef ICE_MODEL
      integer, intent(in) :: Iout, Iuout, Ieout
#  endif
      integer, intent(in) :: Nout
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
# ifdef BIO_GOANPZ
      real(r8), intent(in) :: st(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
#   ifdef LMD_SKPP
      real(r8), intent(in) :: hsbl(LBi:,LBj:)
#   endif
#   ifdef LMD_BKPP
      real(r8), intent(in) :: hbbl(LBi:,LBj:)
#   endif
#   ifdef ICE_MODEL
      real(r8), intent(in) :: ui(LBi:,LBj:,:)
      real(r8), intent(in) :: vi(LBi:,LBj:,:)
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(in) :: hsn(LBi:,LBj:,:)
      real(r8), intent(in) :: tis(LBi:,LBj:)
      real(r8), intent(in) :: ti(LBi:,LBj:,:)
      real(r8), intent(in) :: sfwat(LBi:,LBj:,:)
      real(r8), intent(in) :: ageice(LBi:,LBj:,:)
      real(r8), intent(in) :: io_mflux(LBi:,LBj:)
      real(r8), intent(in) :: sig11(LBi:,LBj:,:)
      real(r8), intent(in) :: sig12(LBi:,LBj:,:)
      real(r8), intent(in) :: sig22(LBi:,LBj:,:)
      real(r8), intent(in) :: t0mk(LBi:,LBj:)
      real(r8), intent(in) :: s0mk(LBi:,LBj:)
      real(r8), intent(in) :: utau_iw(LBi:,LBj:)
      real(r8), intent(in) :: chu_iw(LBi:,LBj:)
#   endif
#   ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
#    ifdef BULK_FLUXES
      real(r8), intent(in) :: lhflx(LBi:,LBj:)
      real(r8), intent(in) :: shflx(LBi:,LBj:)
      real(r8), intent(in) :: lrflx(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#     ifdef EMINUSP
      real(r8), intent(in) :: evap(LBi:,LBj:)
      real(r8), intent(in) :: rain(LBi:,LBj:)
#     endif
#    endif
#    ifdef SHORTWAVE
      real(r8), intent(in) :: srflx(LBi:,LBj:)
#    endif
#   endif
#  endif
#  ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: avgu3d(LBi:,LBj:)
      real(r8), intent(inout) :: avgv3d(LBi:,LBj:)
      real(r8), intent(inout) :: avgt(LBi:,LBj:,:)
# ifdef BIO_GOANPZ
      real(r8), intent(inout) :: avgst(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: avgrho(LBi:,LBj:)
#   ifdef LMD_SKPP
      real(r8), intent(inout) :: avghsbl(LBi:,LBj:)
#   endif
#   ifdef LMD_BKPP
      real(r8), intent(inout) :: avghbbl(LBi:,LBj:)
#   endif
#   ifdef ICE_MODEL
      real(r8), intent(inout) :: avguice(LBi:,LBj:)
      real(r8), intent(inout) :: avgvice(LBi:,LBj:)
      real(r8), intent(inout) :: avgaice(LBi:,LBj:)
      real(r8), intent(inout) :: avghice(LBi:,LBj:)
      real(r8), intent(inout) :: avghsno(LBi:,LBj:)
      real(r8), intent(inout) :: avgtice(LBi:,LBj:)
      real(r8), intent(inout) :: avgtimid(LBi:,LBj:)
      real(r8), intent(inout) :: avgsfwat(LBi:,LBj:)
      real(r8), intent(inout) :: avgageice(LBi:,LBj:)
      real(r8), intent(inout) :: avgiomflx(LBi:,LBj:)
      real(r8), intent(inout) :: avgsig11(LBi:,LBj:)
      real(r8), intent(inout) :: avgsig12(LBi:,LBj:)
      real(r8), intent(inout) :: avgsig22(LBi:,LBj:)
      real(r8), intent(inout) :: avgT0mk(LBi:,LBj:)
      real(r8), intent(inout) :: avgS0mk(LBi:,LBj:)
      real(r8), intent(inout) :: avgutau_iw(LBi:,LBj:)
      real(r8), intent(inout) :: avgchu_iw(LBi:,LBj:)
#   endif
#   ifdef AVERAGES_FLUXES
      real(r8), intent(inout) :: avgstf(LBi:,LBj:)
      real(r8), intent(inout) :: avgswf(LBi:,LBj:)
#    ifdef BULK_FLUXES
      real(r8), intent(inout) :: avglhf(LBi:,LBj:)
      real(r8), intent(inout) :: avgshf(LBi:,LBj:)
      real(r8), intent(inout) :: avglrf(LBi:,LBj:)
      real(r8), intent(inout) :: avguwind(LBi:,LBj:)
      real(r8), intent(inout) :: avgvwind(LBi:,LBj:)
#     ifdef EMINUSP
      real(r8), intent(inout) :: avgevap(LBi:,LBj:)
      real(r8), intent(inout) :: avgrain(LBi:,LBj:)
#     endif
#    endif
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: avgsrf(LBi:,LBj:)
#    endif
#   endif
#  endif
#  ifdef AVERAGES_FLUXES
      real(r8), intent(inout) :: avgsus(LBi:,LBj:)
      real(r8), intent(inout) :: avgsvs(LBi:,LBj:)
      real(r8), intent(inout) :: avgbus(LBi:,LBj:)
      real(r8), intent(inout) :: avgbvs(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: avgu2d(LBi:,LBj:)
      real(r8), intent(inout) :: avgv2d(LBi:,LBj:)
      real(r8), intent(inout) :: avgzeta(LBi:,LBj:)

# else

#  ifdef SOLVE3D
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   ifdef BIO_GOANPZ
      real(r8), intent(in) :: st(LBi:UBi,LBj:UBj,N(ng),3,NTS(ng))
#   endif
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
#   ifdef LMD_SKPP
      real(r8), intent(in) :: hsbl(LBi:UBi,LBj:UBj)
#   endif
#   ifdef LMD_BKPP
      real(r8), intent(in) :: hbbl(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ICE_MODEL
      real(r8), intent(in) :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: tis(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: sfwat(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ageice(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: io_mflux(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sig11(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: sig12(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: sig22(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: t0mk(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s0mk(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: utau_iw(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: chu_iw(LBi:UBi,LBj:UBj)
#   endif
#   ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
#    ifdef BULK_FLUXES
      real(r8), intent(in) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#     ifdef EMINUSP
      real(r8), intent(in) :: evap(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)
#     endif
#    endif
#    ifdef SHORTWAVE
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
#    endif
#   endif
#  endif
#  ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: avgu3d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgv3d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgt(LBi:UBi,LBj:UBj,NT(ng))
#  ifdef BIO_GOANPZ
      real(r8), intent(inout) :: avgst(LBi:UBi,LBJ:UBj,NTS(ng))
#  endif
      real(r8), intent(inout) :: avgrho(LBi:UBi,LBj:UBj)
#   ifdef LMD_SKPP
      real(r8), intent(inout) :: avghsbl(LBi:UBi,LBj:UBj)
#   endif
#   ifdef LMD_BKPP
      real(r8), intent(inout) :: avghbbl(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ICE_MODEL
      real(r8), intent(inout) :: avguice(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgvice(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgaice(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avghice(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avghsno(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgtice(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgtimid(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgsfwat(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgageice(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgiomflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgsig11(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgsig12(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgsig22(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgT0mk(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgS0mk(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgutau_iw(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgchu_iw(LBi:UBi,LBj:UBj)
#   endif
#   ifdef AVERAGES_FLUXES
      real(r8), intent(inout) :: avgstf(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgswf(LBi:UBi,LBj:UBj)
#    ifdef BULK_FLUXES
      real(r8), intent(inout) :: avglhf(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgshf(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avglrf(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avguwind(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgvwind(LBi:UBi,LBj:UBj)
#     ifdef EMINUSP
      real(r8), intent(inout) :: avgevap(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgrain(LBi:UBi,LBj:UBj)
#     endif
#    endif
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: avgsrf(LBi:UBi,LBj:UBj)
#    endif
#   endif
#  endif
#  ifdef AVERAGES_FLUXES
      real(r8), intent(inout) :: avgsus(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgsvs(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgbus(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgbvs(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: avgu2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgv2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgzeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j

      real(r8) :: fac, fac1

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG2(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG2(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG2(ng)).eq.1)).or.                          &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Initialize 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgzeta(i,j)=zeta(i,j,Kout)
            avgu2d (i,j)=ubar(i,j,Kout)
            avgv2d (i,j)=vbar(i,j,Kout)
# ifdef SOLVE3D
#  ifdef LMD_SKPP
            avghsbl(i,j)=hsbl(i,j)
#  endif
#  ifdef LMD_BKPP
            avghbbl(i,j)=hbbl(i,j)
#  endif
#  ifdef ICE_MODEL
            avguice(i,j)=ui(i,j,Iuout)
            avgvice(i,j)=vi(i,j,Iuout)
            avgaice(i,j)=ai(i,j,Iout)
            avghice(i,j)=hi(i,j,Iout)
            avghsno(i,j)=hsn(i,j,Iout)
            avgtice(i,j)=tis(i,j)
            avgtimid(i,j)=ti(i,j,Iout)
            avgsfwat(i,j)=sfwat(i,j,Iout)
            avgageice(i,j)=ageice(i,j,Iout)
            avgiomflx(i,j)=io_mflux(i,j)
            avgsig11(i,j)=sig11(i,j,Ieout)
            avgsig12(i,j)=sig12(i,j,Ieout)
            avgsig22(i,j)=sig22(i,j,Ieout)
            avgT0mk(i,j)=t0mk(i,j)
            avgS0mk(i,j)=s0mk(i,j)
            avgutau_iw(i,j)=utau_iw(i,j)
            avgchu_iw(i,j)=chu_iw(i,j)
#  endif
#  ifdef AVERAGES_FLUXES
            avgstf(i,j)=stflx(i,j,itemp)
            avgswf(i,j)=stflx(i,j,isalt)
#   ifdef BULK_FLUXES
            avglhf(i,j)=lhflx(i,j)
            avgshf(i,j)=shflx(i,j)
            avglrf(i,j)=lrflx(i,j)
            avguwind(i,j)=Uwind(i,j)
            avgvwind(i,j)=Vwind(i,j)
#    ifdef EMINUSP
            avgevap(i,j)=evap(i,j)
            avgrain(i,j)=rain(i,j)
#    endif
#   endif
#   ifdef SHORTWAVE
            avgsrf(i,j)=srflx(i,j)
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
            avgsus(i,j)=sustr(i,j)
            avgsvs(i,j)=svstr(i,j)
            avgbus(i,j)=bustr(i,j)
            avgbvs(i,j)=bvstr(i,j)
# endif
          END DO
        END DO

# ifdef SOLVE3D
!
!  Initialize fields associated with 3D horizontal momentum.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgu3d(i,j)=u(i,j,n(ng),Nout)
            avgv3d(i,j)=v(i,j,n(ng),Nout)
            avgrho(i,j)=rho(i,j,n(ng))
          END DO
        END DO
!
!  Initialized fields associated with tracers.
!
        DO itrc=1,NT(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgt(i,j,itrc)=t(i,j,N(ng),Nout,itrc)
            END DO
          END DO
        END DO
#  ifdef BIO_GOANPZ
        DO itrc=1,NTS(ng)
	    DO j=JstrR,JendR
	      DO i=IstrR,IendR
	        avgst(i,j,itrc)=st(i,j,N(ng),Nout,itrc)
	    END DO
	  END DO
	END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG2(ng)) THEN
!
!  Accumulate 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgzeta(i,j)=avgzeta(i,j)+zeta(i,j,Kout)
            avgu2d (i,j)=avgu2d (i,j)+ubar(i,j,Kout)
            avgv2d (i,j)=avgv2d (i,j)+vbar(i,j,Kout)
# ifdef SOLVE3D
#  ifdef LMD_SKPP
            avghsbl(i,j)=avghsbl(i,j)+hsbl(i,j)
#  endif
#  ifdef LMD_BKPP
            avghbbl(i,j)=avghbbl(i,j)+hbbl(i,j)
#  endif
#  ifdef ICE_MODEL
            avguice(i,j)=avguice(i,j)+ui(i,j,Iuout)
            avgvice(i,j)=avgvice(i,j)+vi(i,j,Iuout)
            avgaice(i,j)=avgaice(i,j)+ai(i,j,Iout)
            avghice(i,j)=avghice(i,j)+hi(i,j,Iout)
            avghsno(i,j)=avghsno(i,j)+hsn(i,j,Iout)
            avgtice(i,j)=avgtice(i,j)+tis(i,j)
            avgtimid(i,j)=avgtimid(i,j)+ti(i,j,Iout)
            avgsfwat(i,j)=avgsfwat(i,j)+sfwat(i,j,Iout)
            avgageice(i,j)=avgageice(i,j)+ageice(i,j,Iout)
            avgiomflx(i,j)=avgiomflx(i,j)+io_mflux(i,j)
            avgsig11(i,j)=avgsig11(i,j)+sig11(i,j,Ieout)
            avgsig12(i,j)=avgsig12(i,j)+sig12(i,j,Ieout)
            avgsig22(i,j)=avgsig22(i,j)+sig22(i,j,Ieout)
            avgT0mk(i,j)=avgT0mk(i,j)+t0mk(i,j)
            avgS0mk(i,j)=avgS0mk(i,j)+s0mk(i,j)
            avgutau_iw(i,j)=avgutau_iw(i,j)+utau_iw(i,j)
            avgchu_iw(i,j)=avgchu_iw(i,j)+chu_iw(i,j)
#  endif
#  ifdef AVERAGES_FLUXES
            avgstf(i,j)=avgstf(i,j)+stflx(i,j,itemp)
            avgswf(i,j)=avgswf(i,j)+stflx(i,j,isalt)
#   ifdef BULK_FLUXES
            avglhf(i,j)=avglhf(i,j)+lhflx(i,j)
            avgshf(i,j)=avgshf(i,j)+shflx(i,j)
            avglrf(i,j)=avglrf(i,j)+lrflx(i,j)
            avguwind(i,j)=avguwind(i,j)+Uwind(i,j)
            avgvwind(i,j)=avgvwind(i,j)+Vwind(i,j)
#    ifdef EMINUSP
            avgevap(i,j)=avgevap(i,j)+evap(i,j)
            avgrain(i,j)=avgrain(i,j)+rain(i,j)
#    endif
#   endif
#   ifdef SHORTWAVE
            avgsrf(i,j)=avgsrf(i,j)+srflx(i,j)
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
            avgsus(i,j)=avgsus(i,j)+sustr(i,j)
            avgsvs(i,j)=avgsvs(i,j)+svstr(i,j)
            avgbus(i,j)=avgbus(i,j)+bustr(i,j)
            avgbvs(i,j)=avgbvs(i,j)+bvstr(i,j)
# endif
          END DO
        END DO

# ifdef SOLVE3D
!
!  Accumulate fields associated with 3D horizontal momentum.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgu3d(i,j)=avgu3d(i,j)+u(i,j,N(ng),Nout)
            avgv3d(i,j)=avgv3d(i,j)+v(i,j,N(ng),Nout)
            avgrho(i,j)=avgrho(i,j)+rho(i,j,N(ng))
          END DO
        END DO
!
!  Accumulate fields associated with tracers.
!
        DO itrc=1,NT(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgt(i,j,itrc)=avgt(i,j,itrc)+t(i,j,N(ng),Nout,itrc)
            END DO
          END DO
        END DO
#  ifdef BIO_GOANPZ
!--------------------------------------------
! values are accumulated in the goanpz.h file
! so do not sum them here
!--------------------------------------------
         DO itrc=1,NTS(ng)
	   DO j=JstrR,JendR
	     DO i=IstrR,IendR
	       avgst(i,j,itrc)=st(i,j,N(ng),Nout,itrc)
	     END DO
	   END DO
	 END DO
#  endif
!
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG2(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG2(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        fac=1.0_r8/REAL(nAVG2(ng),r8)
        IF (SOUTH_WEST_TEST) THEN
          AVG2time(ng)=AVG2time(ng)+REAL(nAVG2(ng),r8)*dt(ng)
        END IF
!
!  Process 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgzeta(i,j)=fac*avgzeta(i,j)
            avgu2d (i,j)=fac*avgu2d (i,j)
            avgv2d (i,j)=fac*avgv2d (i,j)
# ifdef SOLVE3D
#  ifdef LMD_SKPP
            avghsbl(i,j)=fac*avghsbl(i,j)
#  endif
#  ifdef LMD_BKPP
            avghbbl(i,j)=fac*avghbbl(i,j)
#  endif
#  ifdef ICE_MODEL
            avguice(i,j)=fac*avguice(i,j)
            avgvice(i,j)=fac*avgvice(i,j)
            avgaice(i,j)=fac*avgaice(i,j)
            avghice(i,j)=fac*avghice(i,j)
            avghsno(i,j)=fac*avghsno(i,j)
            avgtice(i,j)=fac*avgtice(i,j)
            avgtimid(i,j)=fac*avgtimid(i,j)
            avgsfwat(i,j)=fac*avgsfwat(i,j)
            avgageice(i,j)=fac*avgageice(i,j)
            avgiomflx(i,j)=fac*avgiomflx(i,j)
            avgsig11(i,j)=fac*avgsig11(i,j)
            avgsig12(i,j)=fac*avgsig12(i,j)
            avgsig22(i,j)=fac*avgsig22(i,j)
            avgT0mk(i,j)=fac*avgT0mk(i,j)
            avgS0mk(i,j)=fac*avgS0mk(i,j)
            avgutau_iw(i,j)=fac*avgutau_iw(i,j)
            avgchu_iw(i,j)=fac*avgchu_iw(i,j)
#  endif
#  ifdef AVERAGES_FLUXES
            avgstf(i,j)=fac*avgstf(i,j)
            avgswf(i,j)=fac*avgswf(i,j)
#   ifdef BULK_FLUXES
            avglhf(i,j)=fac*avglhf(i,j)
            avgshf(i,j)=fac*avgshf(i,j)
            avglrf(i,j)=fac*avglrf(i,j)
            avguwind(i,j)=fac*avguwind(i,j)
            avgvwind(i,j)=fac*avgvwind(i,j)
#    ifdef EMINUSP
            avgevap(i,j)=fac*avgevap(i,j)
            avgrain(i,j)=fac*avgrain(i,j)
#    endif
#   endif
#   ifdef SHORTWAVE
            avgsrf(i,j)=fac*avgsrf(i,j)
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
            avgsus(i,j)=fac*avgsus(i,j)
            avgsvs(i,j)=fac*avgsvs(i,j)
            avgbus(i,j)=fac*avgbus(i,j)
            avgbvs(i,j)=fac*avgbvs(i,j)
# endif
          END DO
        END DO

# ifdef SOLVE3D
!
!  Process fields associated with 3D horizontal momentum.
!
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgu3d(i,j)=fac*avgu3d(i,j)
              avgv3d(i,j)=fac*avgv3d(i,j)
              avgrho(i,j)=fac*avgrho(i,j)
          END DO
        END DO
!
!  Process fields associated with tracers.
!
        DO itrc=1,NT(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                avgt(i,j,itrc)=fac*avgt(i,j,itrc)
            END DO
          END DO
        END DO
!#  ifdef BIO_GOANPZ
!-----------------------------------------------------------
! Computes average of time series
! For production, totals are needed so this is commented out
!-----------------------------------------------------------
!        DO itrc=1,NTS(ng)
!	    DO j-JstrR,JendR
!	      DO i-IstrR,IendR
!	        anvst(i,j,itrc)=fac*avgt(i,j,itrc)
!	    END DO
!	  END DO
!	END DO
!#  endif
#  endif
!
      END IF

      RETURN
      END SUBROUTINE set_avg2_tile
#endif
      END MODULE set_avg2_mod
