#include "cppdefs.h"
      MODULE step_pred_mod
#if defined NONLINEAR && defined NEMURO_SAN && defined PREDATOR
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group        John M. Klinck   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-steps  simulated  predators  trajectories using a !
!  fourth-order Milne predictor and fourth-order Hamming corrector.    !
!                                                                      !
!  Vertical diffusion is optionally represented by a random walk,      !
!  in which case a forward scheme is used for vertical displacement.   !
!  The probability distribution for the vertical displacement is       !
!  Gaussian and includes a correction for the vertical gradient in     !
!  diffusion coefficient                                               !
!                                                                      !
! Reference:                                                           !
!                                                                      !
!  Hunter, J.R, P.D. Craig, and H.E. Philips, 1993: On the use of      !
!    random walk models with spatially variable diffusivity,           !
!    Journal of Computational Physics, 106, 366-376.                   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: step_pred

      CONTAINS
!
!***********************************************************************
      SUBROUTINE step_pred (ng, tile, Lstr, Lend)
!***********************************************************************
!
      USE mod_param
      USE mod_pred
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lstr, Lend
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif
      CALL step_pred_tile (ng, tile, Lstr, Lend,                        &
     &                       knew(ng), nnew(ng), nfm3(ng), nfm2(ng),    &
     &                       nfm1(ng), nf(ng), nfp1(ng),                &
     &                       PREDS(ng) % bounded,                       &
     &                       PREDS(ng) % Ftype,                         &
     &                       PREDS(ng) % Tinfo,                         &
     &                       PREDS(ng) % Fz0,                           &
     &                       PREDS(ng) % rwalk,                         &
     &                       PREDS(ng) % bioenergy,                     &
     &                       PREDS(ng) % species,                       &
     &                       PREDS(ng) % alive,                         &
     &                       PREDS(ng) % track)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE step_pred
!
!***********************************************************************
      SUBROUTINE step_pred_tile (ng, tile, Lstr, Lend,                  &
     &                             knew, nnew,                          &
     &                             nfm3, nfm2, nfm1, nf, nfp1,          &
     &                             bounded, Ftype, Tinfo, Fz0,          &
     &                             rwalk, bioenergy, species,           &
     &                             alive, track)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
# ifdef FLOAT_VWALK
      USE mod_mixing
# endif
      USE mod_ncparam
      USE mod_ocean
      USE mod_pred
      USE mod_fish
      USE mod_scalars
      USE mod_biology
      USE pred_swim_mod
      USE interp_floats_mod
!
# ifdef DISTRIBUTE
      USE distribute_mod
# endif
      USE utility_mod, ONLY : nrng
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lstr, Lend
      integer, intent(in) :: knew, nnew, nfm3, nfm2, nfm1, nf, nfp1
!
# ifdef ASSUMED_SHAPE
      integer, intent(in)     :: Ftype(:)
      real(r8), intent(inout) :: Tinfo(0:,:)
      real(r8), intent(inout) :: Fz0(:)
      real(r8), intent(inout) :: rwalk(:)

      logical, intent(inout)  :: bounded(:)
      real(r8), intent(inout) :: bioenergy(:,:)
      integer, intent(inout)  :: species(:)
      logical, intent(inout)  :: alive(:)
      real(r8), intent(inout) :: track(:,0:,:)
# else
      integer, intent(in)     :: Ftype(Npred(ng))
      real(r8), intent(inout) :: Tinfo(0:izrhs,Npred(ng))
      real(r8), intent(inout) :: Fz0(Npred(ng))
      real(r8), intent(inout) :: rwalk(Npred(ng))

      logical, intent(inout)  :: bounded(Npred(ng))
      real(r8), intent(inout) :: bioenergy(NPredV(ng),Npred(ng))
      integer, intent(inout)  :: species(Npred(ng))
      logical, intent(inout)  :: alive(Npred(ng))
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Npred(ng))
# endif
!
!  Local variable declarations.
!
      logical, parameter :: Gmask = .FALSE.
# ifdef MASKING
      logical, parameter :: Lmask = .TRUE.
# else
      logical, parameter :: Lmask = .FALSE.
# endif
      logical, dimension(Lstr:Lend) :: MyPredThread

      integer :: LBi, UBi, LBj, UBj
      integer :: Ir, Jr, Npts, i, i1, i2, j, j1, j2, itrc, l, k
      integer :: NptsF, NptsL, NptsM
      integer :: ipsp, ipid
      integer :: ct

      real(r8), parameter :: Fspv = 0.0_r8
      integer, parameter :: iFspv = 0
      logical, parameter :: lFspv = .false.

      real(r8) :: cff1, cff2, cff3, cff4, cff5, cff6, cff7, cff8, cff9
      real(r8) :: p1, p2, q1, q2, xrhs, yrhs, zrhs, zfloat
      real(r8) :: tspawn, Fworth0, mg_time, xini, yini, zini
!      real(r8) :: Nymort

# ifdef FLOAT_VWALK
      integer :: ierr
      integer :: iseed = 149876
# endif

      real(r8), dimension(Lstr:Lend) :: nudg

# ifdef DISTRIBUTE
      real(r8) :: Xstr, Xend, Ystr, Yend
      real(r8), dimension(Npred(ng)*NFV(ng)*(NFT+1)) :: Fwrk
      real(r8), dimension(Npred(ng)*NPredV(ng)) :: FwrkF
      real(r8), dimension(Npred(ng)*Nspecies(ng)) :: FwrkM
      logical,  dimension(Npred(ng)) :: FwrkL
      integer,  dimension(Npred(ng)) :: FwrkI
# endif
!
! Set tile array bounds.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
! In distributed-memory configuration, determine which node bounds the
! current location of the predator. Assign non-bounded predators to the
! master node.
!-----------------------------------------------------------------------
!
! The strategy here is to build a switch that processes only the pred
! contained within the node bounds. The trajectory data for the new
! time-level (nfp1) is initialized to Fspv. These values are used during
! recombining step at the end of the routine.  Since a SUM reduction is
! carried-out, setting Fspv to zero means the pred only contribute in
! their own tile.
!
      Npts=NFV(ng)*(NFT+1)*Npred(ng)
      NptsF=NPredV(ng)*Npred(ng)
      NptsL=Npred(ng)
      NptsM=Nspecies(ng)*Npred(ng)

      Xstr=REAL(BOUNDS(ng)%Istr(MyRank),r8)-0.5_r8
      Xend=REAL(BOUNDS(ng)%Iend(MyRank),r8)+0.5_r8
      Ystr=REAL(BOUNDS(ng)%Jstr(MyRank),r8)-0.5_r8
      Yend=REAL(BOUNDS(ng)%Jend(MyRank),r8)+0.5_r8
      DO l=Lstr,Lend
        MyPredThread(l)=.FALSE.
        IF ((Xstr.le.track(ixgrd,nf,l)).and.                            &
     &      (track(ixgrd,nf,l).lt.Xend).and.                            &
     &      (Ystr.le.track(iygrd,nf,l)).and.                            &
     &      (track(iygrd,nf,l).lt.Yend)) THEN
          MyPredThread(l)=.TRUE.
        ELSE IF (Master.and.(.not.bounded(l))) THEN
          MyPredThread(l)=.TRUE.
        ELSE
          DO j=0,NFT
            DO i=1,NFV(ng)
              track(i,j,l)=Fspv
            END DO
          END DO
          DO i=1,NPredV(ng)
            bioenergy(i,l)=Fspv
          END DO
          alive(l)=lFspv
        END IF
      END DO
# else
      DO l=Lstr,Lend
        MyPredThread(l)=.TRUE.
      END DO
# endif
!
      DO l=Lstr,Lend
        nudg(l)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Calculate slopes at new time-step.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, ixrhs, Npred(ng),           &
     &                    isUvel, -u3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % u(:,:,:,nnew),                    &
     &                    MyPredThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, iyrhs, Npred(ng),           &
     &                    isUvel, -v3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % v(:,:,:,nnew),                    &
     &                    MyPredThread, bounded, track)

#  if !defined FLOAT_VWALK
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfp1, izrhs, Npred(ng),           &
     &                    isBw3d, w3dvar, Lmask, spval, nudg,           &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % W,                                &
     &                    MyPredThread, bounded, track)
#  endif
# else
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, ixrhs, Npred(ng),           &
     &                    isUbar, -u2dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % ubar(:,:,knew),                   &
     &                    MyPredThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, iyrhs, Npred(ng),           &
     &                    isVbar, -v2dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % vbar(:,:,knew),                   &
     &                    MyPredThread, bounded, track)
# endif
!
!-----------------------------------------------------------------------
!  Predator swimming behavior
!-----------------------------------------------------------------------
!
! If beginning of predator season, randomly initialize position
!      mg_time=REAL(INT(time(ng)/86400.0_r8/days_year))
!      mg_time=time(ng)/86400.0_r8-days_year*mg_time
!      DO l=Lstr,Lend
!        IF (MyPredThread(l).and.bounded(l)) THEN
!          ipsp=idpred(species(l))
!          IF ((mg_time.gt.(Pmgstr(ipsp,ng)-0.5_r8*dt(ng))).and.         &
!     &        (mg_time.lt.(Pmgstr(ipsp,ng)+0.5_r8*dt(ng)))) THEN
!            DO j=0,NFT
!              track(ixgrd,j,l)=Tinfo(ixgrd,l)
!              track(iygrd,j,l)=Tinfo(iygrd,l)
!              track(ixrhs,j,l)=0.0_r8
!              track(iyrhs,j,l)=0.0_r8
!# ifdef SOLVE3D
!              track(izgrd,j,l)=Tinfo(izgrd,l)
!              track(izrhs,j,l)=0.0_r8
!# endif
!            END DO
!          END IF
!        END IF
!      END DO
!
      CALL pred_swim (ng, tile, LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                    nfm1, nf, nfp1, nnew,                         &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
# ifdef SOLVE3D
     &                    GRID(ng) % Hz,                                &
# endif
     &                    MyPredThread, bounded, rwalk, track,          &
     &                    bioenergy, alive, species,                    &
     &                    OCEAN(ng) % pred_count,                       &
     &                    OCEAN(ng) % pred_list,                        &
     &                    PREDS(ng) % prednodes,                        &
     &                    OCEAN(ng) % fish_count,                       &
     &                    OCEAN(ng) % fish_list,                        &
     &                    FISHES(ng) % fishnodes,                       &
     &                    FISHES(ng) % bioenergy,                       &
     &                    FISHES(ng) % species,                         &
     &                    FISHES(ng) % lifestage)
!
!-----------------------------------------------------------------------
!  Determine pred status.
!-----------------------------------------------------------------------
!
# ifdef EW_PERIODIC
      cff1=REAL(Lm(ng),r8)
      DO l=Lstr,Lend
        IF (MyPredThread(l).and.bounded(l)) THEN
          IF (track(ixgrd,nfp1,l).ge.REAL(Lm(ng)+1,r8)-0.5_r8) THEN
            track(ixgrd,nfp1,l)=track(ixgrd,nfp1,l)-cff1
            track(ixgrd,nf  ,l)=track(ixgrd,nf  ,l)-cff1
            track(ixgrd,nfm1,l)=track(ixgrd,nfm1,l)-cff1
            track(ixgrd,nfm2,l)=track(ixgrd,nfm2,l)-cff1
            track(ixgrd,nfm3,l)=track(ixgrd,nfm3,l)-cff1
          ELSE IF (track(ixgrd,nfp1,l).lt.0.5_r8) THEN
            track(ixgrd,nfp1,l)=cff1+track(ixgrd,nfp1,l)
            track(ixgrd,nf  ,l)=cff1+track(ixgrd,nf  ,l)
            track(ixgrd,nfm1,l)=cff1+track(ixgrd,nfm1,l)
            track(ixgrd,nfm2,l)=cff1+track(ixgrd,nfm2,l)
            track(ixgrd,nfm3,l)=cff1+track(ixgrd,nfm3,l)
          END IF
        END IF
      END DO
#  ifdef DISTRIBUTE
      IF (NtileI(ng).gt.1) THEN
        Fwrk=RESHAPE(track,(/Npts/))
        CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
        track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Npred(ng)/))
        DO l=Lstr,Lend
          IF ((Xstr.le.track(ixgrd,nfp1,l)).and.                        &
     &        (track(ixgrd,nfp1,l).lt.Xend).and.                        &
     &        (Ystr.le.track(iygrd,nfp1,l)).and.                        &
     &        (track(iygrd,nfp1,l).lt.Yend)) THEN
            MyPredThread(l)=.TRUE.
          ELSE IF (Master.and.(.not.bounded(l))) THEN
            MyPredThread(l)=.TRUE.
          ELSE
            MyPredThread(l)=.FALSE.
            DO j=0,NFT
              DO i=1,NFV(ng)
                track(i,j,l)=Fspv
              END DO
            END DO
            DO i=1,NPredV(ng)
              bioenergy(i,l)=Fspv
            END DO
            alive(l)=lFspv
          END IF
        END DO
      END IF
#  endif
# else
      DO l=Lstr,Lend
        IF (MyPredThread(l).and.bounded(l)) THEN
! JF: Relocation at starting position for non-periodic case
!          IF ((track(ixgrd,nfp1,l).ge.REAL(Lm(ng)+1,r8)-0.5_r8).or.     &
!     &        (track(ixgrd,nfp1,l).lt.0.5_r8)) THEN
!            DO j=0,NFT
!              track(ixgrd,j,l)=Tinfo(ixgrd,l) 
!              track(iygrd,j,l)=Tinfo(iygrd,l) 
!              track(izgrd,j,l)=Tinfo(izgrd,l) 
!              track(ixrhs,j,l)=0.0_r8
!              track(iyrhs,j,l)=0.0_r8
!              track(izrhs,j,l)=0.0_r8
!            END DO
!          END IF
! JF: Crude reflective boundary condition for non-periodic case
          IF (track(ixgrd,nfp1,l).ge.REAL(Lm(ng)+1,r8)-1.5_r8) THEN
            DO j=0,NFT
              track(ixgrd,j,l)=REAL(Lm(ng)+1-2) 
              track(ixrhs,j,l)=-track(ixrhs,nfp1,l)
            END DO
          END IF
          IF (track(ixgrd,nfp1,l).lt.1.5_r8) THEN
            DO j=0,NFT
              track(ixgrd,j,l)=2.0_r8 
              track(ixrhs,j,l)=-track(ixrhs,nfp1,l)
            END DO
          END IF
! Original ROMS code
!          IF ((track(ixgrd,nfp1,l).ge.REAL(Lm(ng)+1,r8)-0.5_r8).or.     &
!     &        (track(ixgrd,nfp1,l).lt.0.5_r8)) THEN
!            bounded(l)=.FALSE.
!          END IF
        END IF
      END DO
# endif
# ifdef NS_PERIODIC
      cff1=REAL(Mm(ng),r8)
      DO l=Lstr,Lend
        IF (MyPredThread(l).and.bounded(l)) THEN
          IF (track(iygrd,nfp1,l).ge.REAL(Mm(ng)+1,r8)-0.5_r8) THEN
            track(iygrd,nfp1,l)=track(iygrd,nfp1,l)-cff1
            track(iygrd,nf  ,l)=track(iygrd,nf  ,l)-cff1
            track(iygrd,nfm1,l)=track(iygrd,nfm1,l)-cff1
            track(iygrd,nfm2,l)=track(iygrd,nfm2,l)-cff1
            track(iygrd,nfm3,l)=track(iygrd,nfm3,l)-cff1
          ELSE IF (track(iygrd,nfp1,l).lt.0.5_r8) THEN
            track(iygrd,nfp1,l)=cff1+track(iygrd,nfp1,l)
            track(iygrd,nf  ,l)=cff1+track(iygrd,nf  ,l)
            track(iygrd,nfm1,l)=cff1+track(iygrd,nfm1,l)
            track(iygrd,nfm2,l)=cff1+track(iygrd,nfm2,l)
            track(iygrd,nfm3,l)=cff1+track(iygrd,nfm3,l)
          END IF
        END IF
      END DO
#  ifdef DISTRIBUTE
      IF (NtileJ(ng).gt.1) THEN
        Fwrk=RESHAPE(track,(/Npts/))
        CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
        track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Npred(ng)/))
        DO l=Lstr,Lend
          IF ((Xstr.le.track(ixgrd,nfp1,l)).and.                        &
     &        (track(ixgrd,nfp1,l).lt.Xend).and.                        &
     &        (Ystr.le.track(iygrd,nfp1,l)).and.                        &
     &        (track(iygrd,nfp1,l).lt.Yend)) THEN
            MyPredThread(l)=.TRUE.
          ELSE IF (Master.and.(.not.bounded(l))) THEN
            MyPredThread(l)=.TRUE.
          ELSE
            MyPredThread(l)=.FALSE.
            DO j=0,NFT
              DO i=1,NFV(ng)
                track(i,j,l)=Fspv
              END DO
            END DO
            DO i=1,NPredV(ng)
                bioenergy(i,l)=Fspv
            END DO
            alive(l)=lFspv
          END IF
        END DO
      END IF
#  endif
# else
      DO l=Lstr,Lend
        IF (MyPredThread(l).and.bounded(l)) THEN
! JF: Relocation at starting position for non-periodic case
!          IF ((track(iygrd,nfp1,l).ge.REAL(Mm(ng)+1,r8)-0.5_r8).or.     &
!     &        (track(iygrd,nfp1,l).lt.0.5_r8)) THEN
!            DO j=0,NFT
!              track(ixgrd,j,l)=Tinfo(ixgrd,l) 
!              track(iygrd,j,l)=Tinfo(iygrd,l) 
!              track(izgrd,j,l)=Tinfo(izgrd,l) 
!              track(ixrhs,j,l)=0.0_r8
!              track(iyrhs,j,l)=0.0_r8
!              track(izrhs,j,l)=0.0_r8
!            END DO
!          END IF
! JF: Crude reflective boundary condition for non-periodic case
          IF (track(iygrd,nfp1,l).ge.REAL(Mm(ng)+1,r8)-1.5_r8) THEN
            DO j=0,NFT
              track(iygrd,j,l)=REAL(Mm(ng)+1-2) 
              track(iyrhs,j,l)=-track(iyrhs,nfp1,l)
            END DO
          END IF
          IF (track(iygrd,nfp1,l).lt.1.5_r8) THEN
            DO j=0,NFT
              track(iygrd,j,l)=2.0_r8 
              track(iyrhs,j,l)=-track(iyrhs,nfp1,l)
            END DO
          END IF
! Original ROMS code
!          IF ((track(iygrd,nfp1,l).ge.REAL(Mm(ng)+1,r8)-0.5_r8).or.     &
!     &        (track(iygrd,nfp1,l).lt.0.5_r8)) THEN
!            bounded(l)=.FALSE.
!          END IF
        END IF
      END DO
# endif
# ifdef SOLVE3D
!
!  Reflect pred at surface or bottom.
!
      DO l=Lstr,Lend
        IF (MyPredThread(l).and.bounded(l)) THEN
          IF (track(izgrd,nfp1,l).gt.REAL(N(ng),r8))                    &
     &      track(izgrd,nfp1,l)=2.0_r8*REAL(N(ng),r8)-                  &
     &                          track(izgrd,nfp1,l)
          IF (track(izgrd,nfp1,l).lt.0.0_r8)                            &
     &      track(izgrd,nfp1,l)=-track(izgrd,nfp1,l)
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Calculate slopes with corrected locations.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, ixrhs, Npred(ng),           &
     &                    isUvel, -u3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % u(:,:,:,nnew),                    &
     &                    MyPredThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, iyrhs, Npred(ng),           &
     &                    isVvel, -v3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % v(:,:,:,nnew),                    &
     &                    MyPredThread, bounded, track)

#  if !defined FLOAT_VWALK
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfp1, izrhs, Npred(ng),           &
     &                    isBw3d, -w3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % W,                                &
     &                    MyPredThread, bounded, track)
#  endif
# else
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, ixrhs, Npred(ng),           &
     &                    isUbar, -u2dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % ubar(:,:,knew),                   &
     &                    MyPredThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, iyrhs, Npred(ng),           &
     &                    isVbar, -v2dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % vbar(:,:,knew),                   &
     &                    MyPredThread, bounded, track)
# endif
!
!  If newly released pred, initialize slopes at all time levels.
!
!      DO l=Lstr,Lend
!        IF (MyPredThread(l).and.bounded(l).and.                             &
!     &      (time(ng)-dt(ng).le.Tinfo(itstr,l).and.                     &
!     &       time(ng)+dt(ng).gt.Tinfo(itstr,l))) THEN
!          xrhs=track(ixrhs,nfp1,l)
!          yrhs=track(iyrhs,nfp1,l)
!# ifdef SOLVE3D
!          zrhs=track(izrhs,nfp1,l)
!# endif
!          DO i=0,NFT
!            track(ixrhs,i,l)=xrhs
!            track(iyrhs,i,l)=yrhs
!# ifdef SOLVE3D
!            track(izrhs,i,l)=zrhs
!# endif
!          END DO
!        END IF
!      END DO
!
!-----------------------------------------------------------------------
!  Interpolate various output variables at the corrected locations.
!-----------------------------------------------------------------------
!
      IF (spherical) THEN
        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflon, Npred(ng),         &
     &                      isBr2d, r2dvar, Gmask, spval, nudg,         &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % lonr,                            &
     &                      MyPredThread, bounded, track)

        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflat, Npred(ng),         &
     &                      isBr2d, r2dvar, Gmask, spval, nudg,         &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % latr,                            &
     &                      MyPredThread, bounded, track)
      ELSE
        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflon, Npred(ng),         &
     &                      isBr2d, r2dvar, Gmask, spval, nudg,         &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % xr,                              &
     &                      MyPredThread, bounded, track)

        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflat, Npred(ng),         &
     &                      isBr2d, r2dvar, Gmask, spval, nudg,         &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % yr,                              &
     &                      MyPredThread, bounded, track)
      END IF
# ifdef SOLVE3D

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfp1, idpth, Npred(ng),           &
     &                    isBw3d, w3dvar, Lmask, spval, nudg,           &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    GRID(ng) % z_w,                               &
     &                    MyPredThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, ifden, Npred(ng),           &
     &                    isBr3d, r3dvar, Lmask, spval, nudg,           &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % rho,                              &
     &                    MyPredThread, bounded, track)

      DO itrc=1,NT(ng)
        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                      Lstr, Lend, nfp1, itrc+NFV(ng)-NT(ng),      &
     &                      Npred(ng), isTvar(itrc),                    &
     &                      r3dvar, Lmask, spval, nudg,                 &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % Hz,                              &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      OCEAN(ng) % t(:,:,:,nnew,itrc),             &
     &                      MyPredThread, bounded, track)
      END DO
!
# endif
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect pred on all nodes.
!-----------------------------------------------------------------------
!
      !Fwrk=RESHAPE(track,(/Npts/))
      !CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
      !track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Npred(ng)/))

      ct=1
      DO l=1,Npred(ng)
        DO j=0,NFT
          DO i=1,NFV(ng)
               Fwrk(ct) = track(i,j,l)
               ct=ct+1
          ENDDO
        ENDDO
      ENDDO

      CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)

      ct=1
      DO l=1,Npred(ng)
        DO j=0,NFT
          DO i=1,NFV(ng)
               track(i,j,l) = Fwrk(ct)
               ct=ct+1
          ENDDO
        ENDDO
      ENDDO

      !FwrkF=RESHAPE(bioenergy,(/NptsF/))
      !CALL mp_collect (ng, iNLM, NptsF, Fspv, FwrkF)
      !bioenergy=RESHAPE(FwrkF,(/NPredV(ng),Npred(ng)/))

      ct=1
      DO l=1,Npred(ng)
        DO j=1,NpredV(ng)
            FwrkF(ct) = bioenergy(j,l)
            ct=ct+1
        ENDDO
      ENDDO

      CALL mp_collect (ng, iNLM, NptsF, Fspv, FwrkF)

      ct=1
      DO l=1,Npred(ng)
        DO j=1,NpredV(ng)
            bioenergy(j,l) = FwrkF(ct)
            ct=ct+1
        ENDDO
      ENDDO

      CALL mp_collect_l (ng, iNLM, Npred(ng), alive)
!
#ifdef FOO
!  Collect the alive status.
      Fwrk=Fspv
      DO l=1,Npred(ng)
        IF (bounded(l)) THEN
          Fwrk(l)=1.0_r8
        END IF
      END DO
      CALL mp_collect (ng, iNLM, Npred(ng), Fspv, Fwrk)
      DO l=1,Npred(ng)
        IF (Fwrk(l).ne.Fspv) THEN
          alive(l)=.TRUE.
        ELSE
          alive(l)=.FALSE.
        END IF
      END DO
#endif
!
!  Collect the bounded status switch.
      Fwrk=Fspv
      DO l=1,Npred(ng)
        IF (bounded(l)) THEN
          Fwrk(l)=1.0_r8
        END IF
      END DO
      CALL mp_collect (ng, iNLM, Npred(ng), Fspv, Fwrk)
      DO l=1,Npred(ng)
        IF (Fwrk(l).ne.Fspv) THEN
          bounded(l)=.TRUE.
        ELSE
          bounded(l)=.FALSE.
        END IF
      END DO
!
# endif
      RETURN
      END SUBROUTINE step_pred_tile
#endif
      END MODULE step_pred_mod
