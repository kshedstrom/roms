#include "cppdefs.h"
      MODULE pred_swim_mod
#if defined NONLINEAR && defined NEMURO_SAN && defined PREDATOR
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group         Mark Hadfield   !
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
!     fishthread   Float parallel thread bounded switch.                 !
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
      SUBROUTINE pred_swim (ng, tile, LBi, UBi, LBj, UBj, LBk, UBk,     &
     &                          itimem1, itime, itimep1, nnew,          &
     &                          pm, pn,                                 &
# ifdef SOLVE3D
     &                          Hz,                                     &
# endif
     &                          fishthread, bounded, rwalk, track,        &
     &                          bioenergy, alive, species,              &
     &                          pred_count, pred_list, prednodes,       &
     &                          fish_count, fish_list, fishnodes,       &
     &                          fish_bioenergy, fish_species,           &
     &                          fish_lifestage)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
      USE mod_fish
      USE mod_grid
      USE mod_parallel
      USE mod_types
      USE mod_parallel
      USE nrutil
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: itimem1, itime, itimep1, nnew
      integer, intent(in) :: species(Npred(ng))
      integer, intent(in) :: fish_species(Nfish(ng))
      integer, intent(in) :: pred_count(LBi:UBi,LBj:UBj)
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)
      integer, intent(in) :: fish_lifestage(Nfish(ng))

      type(prednode), intent(in) :: pred_list(LBi:UBi,LBj:UBj)
      type(prednode), target, intent(in) :: prednodes(Npred(ng))
      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Nfish(ng))

      logical, intent(in) :: fishthread(Npred(ng))
      logical, intent(in) :: bounded(Npred(ng))
      logical, intent(in) :: alive(Npred(ng))

      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
# endif
      real(r8), intent(in) :: bioenergy(NPredV(ng),Npred(ng))
      real(r8), intent(inout) :: rwalk(Npred(ng)*6)
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Npred(ng))
      real(r8), intent(in) :: fish_bioenergy(NFishV(ng),Nfish(ng))
!
!  Local variable declarations.
!
      integer :: Ir, Jr, Kr, IrMax, JrMax
      integer :: i1, i2, j1, j2, i, j, ii, jj, itm1
      integer :: ipred, ipid, ipsp, ifish, ifid, ifsp
      integer :: iopt, iind(3), jind(3), nind

      type(prednode), pointer :: thispred
      type(fishnode), pointer :: thisfish

      real(r8) :: Fweight, Fworth, score, scoreMax
      real(r8) :: d_i, d_j, d_ij, theta_ij, fac, snudg
      real(r8) :: mg_time, xini, yini, zini, dtpred

# include "set_bounds.h"
!
# ifdef DISTRIBUTE
      IF (Master) THEN
        CALL ran1 (rwalk)
      END IF
      CALL mp_bcastf (ng, iNLM, rwalk)
# elif defined _OPENMP
!$OMP SINGLE
      CALL ran1 (rwalk)
!$OMP END SINGLE
# else
!     IF (Lstr.eq.1) THEN
        CALL ran1 (rwalk)
!     END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute swimming velocity based on fitness behavior (Railsback)
!-----------------------------------------------------------------------
!
!      dtpred=dt(ng)
      dtpred=3600.0_r8  ! hourly
!
! mg_time: time in days, modulo days_year
! KLUDGE: Depends on time origin being at the start of a year
      mg_time=REAL(INT(time(ng)/86400.0_r8/days_year))
      mg_time=time(ng)/86400.0_r8-days_year*mg_time
!
!      DO i=LBi+1,UBi-1
!        DO j=LBj+1,UBj-1
      DO i=MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
        DO j=MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
          IF (pred_count(i,j).gt.0) THEN
            thispred => pred_list(i,j) % next
            DO ipred=1,pred_count(i,j)
              ipid = thispred % pred
              ipsp = idpred(species(ipid))
              IF (fishthread(ipid).and.bounded(ipid)                      &
     &            .and.alive(ipid)) THEN
                IF ((mg_time.ge.Pmgstr(ipsp,ng)).and.                   &
     &              (mg_time.le.Pmgend(ipsp,ng))) THEN
!
                  Ir=NINT(track(ixgrd,itime,ipid))
                  Ir=MIN(MAX(Ir,1),Lm(ng))
                  Jr=NINT(track(iygrd,itime,ipid))
                  Jr=MIN(MAX(Jr,1),Mm(ng))

                  scoreMax=-99.0_r8
                  IrMax=Ir
                  JrMax=Jr
!
! Randomization of i indices
!                 iopt=INT(6.0_r8*rwalk(ipid))+1
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
!                 iopt=INT(6.0_r8*rwalk(ipid+Npred(ng)))+1
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
                    score=0.0_r8
                    IF (fish_count(ii,jj).gt.0) THEN
                      thisfish => fish_list(ii,jj) % next
                      DO ifish=1,fish_count(ii,jj)
                        ifid = thisfish % fish
                        ifsp = idfish(fish_species(ifid))
                        IF (fish_lifestage(ifid).ge.if_subadult) THEN
                          Fweight=fish_bioenergy(ifwwt,ifid)
                          Fworth=fish_bioenergy(ifworth,ifid)
                          score=score+Fweight*Fworth*                   &
     &                                Fpref(ifsp,ng)/K_Fish(ifsp,ng)
                        END IF
                        thisfish => thisfish % next
                      END DO
                    END IF
                    IF (score.gt.scoreMax) THEN
                      scoreMax=score
                      IrMax=ii
                      JrMax=jj
                    END IF
                  END DO
!                fac=2.0_r8*(rwalk(ipid+2*Npred(ng))-0.5_r8)
                  CALL ran1 (snudg)
                  fac=2.0_r8*(snudg-0.5_r8)
                  d_i=REAL(IrMax,r8)-track(ixgrd,itime,ipid)
                  d_j=REAL(JrMax,r8)-track(iygrd,itime,ipid)
! Angle with +-0.5 radians randomness
                  theta_ij=atan2(d_j,d_i)+0.5_r8*fac
! Distance with 30% randomness
                  d_ij=((d_i/pm(Ir,Jr))**2+(d_j/pn(Ir,Jr))**2)**0.5_r8
! Pswim is in cm/s, so multiply by 0.01 for m/s)
                  d_ij=MIN(0.01_r8*Pswim(ipsp,ng)*dtpred,d_ij)
                  d_ij=d_ij+0.3_r8*fac*d_ij
! Swim towards best location
                  track(ixgrd,itimep1,ipid)=track(ixgrd,itime,ipid)+    &
     &                                     d_ij*cos(theta_ij)*pm(Ir,Jr)
                  track(iygrd,itimep1,ipid)=track(iygrd,itime,ipid)+    &
     &                                     d_ij*sin(theta_ij)*pn(Ir,Jr)
                ELSE
! Freeze location (outside of migration season)
                  track(ixgrd,itimep1,ipid)=track(ixgrd,itime,ipid)
                  track(iygrd,itimep1,ipid)=track(iygrd,itime,ipid)
                END IF
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
