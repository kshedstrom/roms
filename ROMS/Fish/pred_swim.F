#include "cppdefs.h"
      MODULE pred_swim_mod
#if defined NONLINEAR && defined NEMURO_SAN
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group         Mark Hadfield   !
!    Licensed under a MIT/X style license             John M. Klinck   !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine interpolates requested field at the float trajectory   !
!  locations.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     Lstr       Starting float index to process.                      !
!     Lend       Ending   float index to process.                      !
!     itime      Floats time level to process.                         !
!     ifield     ID of field to compute.                               !
!     gtype      Grid type. If negative, interpolate floats slopes.    !
!     maskit     Should the field be masked? Ignored if Land/Sea       !
!                 masking is not active.                               !
!     nudg       Vertical random walk term to be added to the field.   !
!     pm         Inverse grid spacing (1/m) in the XI-direction.       !
!     pn         Inverse grid spacing (1/m) in the ETA-direction.      !
!     Hz         Vertical thicknesses (m).                             !
!     Amask      Field Land/Sea mask.                                  !
!     A          Field to interpolate from.                            !
!     MyThread   Float parallel thread bounded switch.                 !
!     bounded    Float grid bounded status switch.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     track      Interpolated field: track(ifield,itime,:).            !
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: pred_swim

      CONTAINS
!
!***********************************************************************
      SUBROUTINE pred_swim (ng, LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          itimem1, itime, itimep1, nnew,          &
     &                          pm, pn,                                 &
# ifdef SOLVE3D
     &                          Hz,                                     &
# endif
     &                          MyThread, bounded, track,               &
     &                          bioenergy, alive, species,              &
     &                          pred_count, pred_list, prednodes,       &
     &                          fish_count, fish_list, fishnodes)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
      USE mod_grid
      USE mod_types
      USE utility_mod, ONLY : nrng
      USE nrutil
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: itimem1, itime, itimep1, nnew
      integer, intent(in) :: species(Npred(ng))
      integer, intent(in) :: pred_count(LBi:UBi,LBj:UBj)
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)

      type(prednode), intent(in) :: pred_list(LBi:UBi,LBj:UBj)
      type(prednode), target, intent(in) :: prednodes(Npred(ng))
      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Nfish(ng))

      logical, intent(in) :: MyThread(Npred(ng))
      logical, intent(in) :: bounded(Npred(ng))
      logical, intent(in) :: alive(Npred(ng))

      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
# endif
      real(r8), intent(in) :: bioenergy(NPredV(ng),Npred(ng))
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Npred(ng))
!
!  Local variable declarations.
!
      integer :: Ir, Jr, Kr, IrMax, JrMax
      integer :: i1, i2, j1, j2, i, j, ii, jj, itm1
      integer :: ipred, pid, isp, ifish, fid
      integer :: iopt, iind(3), jind(3), nind

      type(prednode), pointer :: thispred
      type(fishnode), pointer :: thisfish

      real(r8) :: Pweight, score, scoreMax, Uswim
      real(r8) :: Fweight, Fworth, Fbiomass
      real(r8) :: d_i, d_j, d_ij, theta_ij, fac, snudg
!
!-----------------------------------------------------------------------
!  Compute swimming velocity based on Railsback or Humston behavior
!-----------------------------------------------------------------------
!
      DO i=LBi+1,UBi-1
        DO j=LBj+1,UBj-1
          IF (pred_count(i,j).gt.0) THEN
            thispred => pred_list(i,j) % next
            DO ipred=1,pred_count(i,j)
              pid = thispred % pred
              isp = idpred(species(pid))
              IF (MyThread(pid).and.bounded(pid).and.alive(pid)) THEN
!
                Ir=NINT(track(ixgrd,itime,pid))
                Ir=MIN(MAX(Ir,1),Lm(ng))
                Jr=NINT(track(iygrd,itime,pid))
                Jr=MIN(MAX(Jr,1),Mm(ng))
                Pweight=bioenergy(ipwwt,pid)

! Swimming velocity based on fish weight (in cm/s)
                Uswim=a_Swim(isp,ng)*Pweight**b_Swim(isp,ng)
!
! -------- FITNESS-BASED SWIMMING BEHAVIOR (RAILSBACK) --------
!
                scoreMax=-99.0_r8
                IrMax=Ir
                JrMax=Jr
!
! Randomization of i indices
                CALL ran1 (snudg)
                iopt=INT(6.0_r8*snudg)+1
                IF (iopt.eq.1) THEN
                  iind(1)=MAX(Ir-1,1)
                  iind(2)=Ir
                  iind(3)=MIN(Ir+1,Lm(ng))
                ELSE IF (iopt.eq.2) THEN
                  iind(1)=MAX(Ir-1,1)
                  iind(3)=Ir
                  iind(2)=MIN(Ir+1,Lm(ng))
                ELSE IF (iopt.eq.3) THEN
                  iind(2)=MAX(Ir-1,1)
                  iind(1)=Ir
                  iind(3)=MIN(Ir+1,Lm(ng))
                ELSE IF (iopt.eq.4) THEN
                  iind(2)=MAX(Ir-1,1)
                  iind(3)=Ir
                  iind(1)=MIN(Ir+1,Lm(ng))
                ELSE IF (iopt.eq.5) THEN
                  iind(3)=MAX(Ir-1,1)
                  iind(1)=Ir
                  iind(2)=MIN(Ir+1,Lm(ng))
                ELSE IF (iopt.eq.6) THEN
                  iind(3)=MAX(Ir-1,1)
                  iind(2)=Ir
                  iind(1)=MIN(Ir+1,Lm(ng))
                END IF
! Randomization of j indices
                CALL ran1 (snudg)
                iopt=INT(6.0_r8*snudg)+1
                IF (iopt.eq.1) THEN
                  jind(1)=MAX(Jr-1,1)
                  jind(2)=Jr
                  jind(3)=MIN(Jr+1,Mm(ng))
                ELSE IF (iopt.eq.2) THEN
                  jind(1)=MAX(Jr-1,1)
                  jind(3)=Jr
                  jind(2)=MIN(Jr+1,Mm(ng))
                ELSE IF (iopt.eq.3) THEN
                  jind(2)=MAX(Jr-1,1)
                  jind(1)=Jr
                  jind(3)=MIN(Jr+1,Mm(ng))
                ELSE IF (iopt.eq.4) THEN
                  jind(2)=MAX(Jr-1,1)
                  jind(3)=Jr
                  jind(1)=MIN(Jr+1,Mm(ng))
                ELSE IF (iopt.eq.5) THEN
                  jind(3)=MAX(Jr-1,1)
                  jind(1)=Jr
                  jind(2)=MIN(Jr+1,Mm(ng))
                ELSE IF (iopt.eq.6) THEN
                  jind(3)=MAX(Jr-1,1)
                  jind(2)=Jr
                  jind(1)=MIN(Jr+1,Mm(ng))
                END IF
!
                DO nind=1,3
                  ii=iind(nind)
                  jj=jind(nind)
                  Fbiomass=0.0_r8
                  IF (fish_count(i,j).gt.0) THEN
                    thisfish => fish_list(i,j) % next
                    DO ifish=1,fish_count(i,j)
                      fid = thisfish % fish
                      Fweight=bioenergy(ifwwt,pid)
                      Fworth=bioenergy(ifworth,pid)
                      Fbiomass = Fbiomass+Fweight*Fworth
                      thisfish => thisfish % next
                    END DO
                  END IF
                  score=exp(Fbiomass)
                  IF (score.gt.scoreMax) THEN
                    scoreMax=score
                    IrMax=ii
                    JrMax=jj
                  END IF
                END DO
                CALL ran1 (snudg)
                fac=2.0_r8*(snudg-0.5_r8)
                d_i=REAL(IrMax,r8)-track(ixgrd,itime,fid)
                d_j=REAL(JrMax,r8)-track(iygrd,itime,fid)
! Angle with +-0.5 radians randomness
                theta_ij=atan2(d_j,d_i)+0.5_r8*fac
! Distance with 30% randomness
                d_ij=((d_i/pm(Ir,Jr))**2+(d_j/pn(Ir,Jr))**2)**0.5_r8
! Uswim is in cm/s, so multiply by 0.01 for m/s)
                d_ij=MIN(0.01_r8*Uswim*dt(ng),d_ij)
                d_ij=d_ij+0.3_r8*fac*d_ij
! Swim towards best location
                track(ixgrd,itimep1,fid)=track(ixgrd,itime,fid)+      &
     &                                   d_ij*cos(theta_ij)*pm(Ir,Jr)
                track(iygrd,itimep1,fid)=track(iygrd,itime,fid)+      &
     &                                   d_ij*sin(theta_ij)*pn(Ir,Jr)
!
              END IF
              thispred => thispred % next
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE pred_swim
#endif
      END MODULE pred_swim_mod
