#include "cppdefs.h"
      MODULE fish_swim_mod
#if defined NONLINEAR && defined NEMURO_SAN
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
      PUBLIC  :: fish_swim

      CONTAINS
!
!***********************************************************************
      SUBROUTINE fish_swim (ng, tile, LBi, UBi, LBj, UBj, LBk, UBk,     &
     &                          itimem2, itimem1, itime, itimep1, nnew, &
     &                          pm, pn, h,                              &
# ifdef MASKING
     &                          rmask,                                  &
# endif
# ifdef SOLVE3D
     &                          Hz,                                     &
# endif
     &                          fishthread, bounded, rwalk, r2walk,       &
     &                          track, bioenergy, alive,                &
     &                          species, lifestage, swimtype,           &
     &                          fish_count, fish_list, fishnodes)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
      USE mod_grid
      USE mod_types
      USE nrutil
      USE mod_fish
      USE mod_ocean
      USE mod_parallel
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: itimem2, itimem1, itime, itimep1, nnew
      integer, intent(in) :: species(Nfish(ng))
      integer, intent(in) :: lifestage(Nfish(ng))
      integer, intent(in) :: swimtype(2,Nfish(ng))
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)

      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Nfish(ng))

      logical, intent(in) :: fishthread(Nfish(ng))
      logical, intent(in) :: bounded(Nfish(ng))
      logical, intent(in) :: alive(Nfish(ng))

      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
# endif
      real(r8), intent(in) :: bioenergy(NFishV(ng),Nfish(ng))
      real(r8), intent(inout) :: rwalk(Nfish(ng)*3)
      real(r8), intent(inout) :: r2walk(Nfish(ng)*6)
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfish(ng))
!
!  Local variable declarations.
!
      logical :: lflag=.TRUE.

      integer :: ierr
      integer :: iseed = 149876
      integer :: Ir, Jr, Kr, IrMax, JrMax, KrMax, kZmax
      integer :: i1, i2, j1, j2, i, j, ii, jj, kk, itm1
      integer :: ifish, fid, isp, ils
      integer :: iopt, iind(3), jind(3), nind

      type(fishnode), pointer :: thisfish

      real(r8), parameter :: MinVal = 1.0e-6_r8

      real(r8) :: SZooC, LZooC, PZooC, TZooC, ZooMax
      real(r8) :: Fweight, Fgrowth, Flength
      real(r8) :: Fmort
      real(r8) :: score, scoreMax, Uswim
      real(r8) :: d_i, d_j, d_ij, d_k, theta_ij, fac, snudg
      real(r8) :: Fvar, Fvar_opt, Fvar_sigma, zf_lon, zf_lat
      real(r8) :: fTemp, fTemp_opt, fTemp_sig, fTemp_alpha, facT
      real(r8) :: fPval, fPval_opt, fPval_sig, fPval_alpha, facPv
      real(r8) :: cue1, cue1_f, cue1_g, cue1_a
      real(r8) :: cue2, cue2_f, cue2_g, cue2_a
      real(r8) :: dist, epsd, epsd_sigma, ee, epsx, epsy, h1, h2, tau
      real(r8) :: f_xdis, g_epsx, f_ydis, g_epsy, xdis, ydis
      real(r8) :: f_zdis, g_epsz, epsz, zdis, zpos, cff1, cff2, zkdis
      real(r8) :: Bsz, Blz, Bpz, Csz, Clz, Cpz
      real(r8) :: tcheck, mg_time, dtfish

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize some random numbers
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      IF (Master) THEN
#  ifdef CLAMPED
        CALL gasdev_clamp (rwalk, 3.0_r8)
#  else
        DO i=1,SIZE(rwalk)
          IF (rwalk(i) > 3.0_r8) rwalk(i) = 3.0_r8
          IF (rwalk(i) < -3.0_r8) rwalk(i) = -3.0_r8
        END DO
#  endif
        CALL ran1 (r2walk)
      END IF
      CALL mp_bcastf (ng, iNLM, rwalk)
      CALL mp_bcastf (ng, iNLM, r2walk)
# elif defined _OPENMP
!$OMP SINGLE
#  ifdef CLAMPED
      CALL gasdev_clamp (rwalk, 3.0_r8)
#  else
      DO i=1,SIZE(rwalk)
        IF (rwalk(i) > 3.0_r8) rwalk(i) = 3.0_r8
        IF (rwalk(i) < -3.0_r8) rwalk(i) = -3.0_r8
      END DO
#  endif
      CALL ran1 (r2walk)
!$OMP END SINGLE
# else
!     IF (Lstr.eq.1) THEN
#  ifdef CLAMPED
        CALL gasdev_clamp (rwalk, 3.0_r8)
#  else
        DO i=1,SIZE(rwalk)
          IF (rwalk(i) > 3.0_r8) rwalk(i) = 3.0_r8
          IF (rwalk(i) < -3.0_r8) rwalk(i) = -3.0_r8
        END DO
#  endif
        CALL ran1 (r2walk)
!     END IF
# endif
!
!
!-----------------------------------------------------------------------
!  Compute swimming velocity based on Railsback or Humston behavior
!-----------------------------------------------------------------------
!
      dtfish = 3600.0_r8  ! one hour
!
! Note: diurnal vertcial magration for larva could be imposed as:
!       depth = (Fz0/2)*(1+cos(dt*2*pi/86400))
!
!      DO i=LBi+1,UBi-1
!        DO j=LBj+1,UBj-1
      DO i=MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
        DO j=MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              fid = thisfish % fish
              isp = idfish(species(fid))
              ils = lifestage(fid)
              IF (fishthread(fid).and.bounded(fid).and.                   &
     &            alive(fid).and.(lifestage(fid).ge.if_larva)) THEN
!
                Ir=NINT(track(ixgrd,itime,fid))
                Ir=MIN(MAX(Ir,1),Lm(ng))
                Jr=NINT(track(iygrd,itime,fid))
                Jr=MIN(MAX(Jr,1),Mm(ng))
                Kr=NINT(track(izgrd,itime,fid))
                Kr=MIN(MAX(Kr,2),N(ng))
                Flength=bioenergy(iflngth,fid)
                fTemp=track(itemp+NFV(ng)-NT(ng),itime,fid)
! Swimming velocity based on body length per second (in cm/s)
                Uswim=0.1_r8*Fswim(ils,isp,ng)*Flength
!
!
! -------- NO MOVEMENT (ONLY FOR DEBUGGING) --------
                IF ((swimtype(1,fid).eq.0).and.                         &
     &              (swimtype(2,fid).eq.0)) THEN
!
                  track(ixgrd,itimep1,fid)=track(ixgrd,itime,fid)
                  track(iygrd,itimep1,fid)=track(iygrd,itime,fid)
!                  track(izgrd,itimep1,fid)=track(izgrd,itime,fid)
                  track(izgrd,itimep1,fid)=REAL(N(ng),r8)
!
! -------- FITNESS-BASED SWIMMING BEHAVIOR (RAILSBACK) --------
                ELSE IF ((swimtype(1,fid).eq.1).and.                    &
     &                   (swimtype(2,fid).eq.1)) THEN
!
                  scoreMax=-99.0_r8
                  IrMax=Ir
                  JrMax=Jr
                  KrMax=Kr
!
! Randomization of i indices
!                  iopt=INT(6.0_r8*r2walk(fid+3*Nfish(ng)))+1
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
!                  iopt=INT(6.0_r8*r2walk(fid+4*Nfish(ng)))+1
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
                  DO nind=1,3
                    ii=iind(nind)
                    jj=jind(nind)
                    DO kk=N(ng),2,-1
                      CALL fish_growth_ijk (ng, ii, jj, kk, nnew, isp,  &
     &                   ils, Fweight, Flength, fTemp, Uswim, Fgrowth)
! MISSING HERE: CALL FOR MORTALITY SUBROUTINE
                      Fmort=0.0_r8
                      score=exp(Fgrowth)*exp(Fmort)
                      IF (score.gt.scoreMax) THEN
                        scoreMax=score
                        IrMax=ii
                        JrMax=jj
                        KrMax=kk
                      END IF
                    END DO
                  END DO
!                  fac=2.0_r8*(r2walk(fid+5*Nfish(ng))-0.5_r8)
                  CALL ran1 (snudg)
                  fac=2.0_r8*(snudg-0.5_r8)
                  d_i=REAL(IrMax,r8)-track(ixgrd,itime,fid)
                  d_j=REAL(JrMax,r8)-track(iygrd,itime,fid)
! Angle with +-0.5 radians randomness
                  theta_ij=atan2(d_j,d_i)+0.5_r8*fac
! Distance with 30% randomness
                  d_ij=((d_i/pm(Ir,Jr))**2+(d_j/pn(Ir,Jr))**2)**0.5_r8
! Uswim is in cm/s, so multiply by 0.01 for m/s)
                  d_ij=MIN(0.01_r8*Uswim*dtfish,d_ij)
                  d_ij=d_ij+0.3_r8*fac*d_ij
! Swim towards best location
                  IF (lifestage(fid).eq.if_larva) THEN
                    track(ixgrd,itimep1,fid)=track(ixgrd,itime,fid)
                    track(iygrd,itimep1,fid)=track(iygrd,itime,fid)
                  ELSE
                    track(ixgrd,itimep1,fid)=track(ixgrd,itime,fid)+    &
     &                                  d_ij*cos(theta_ij)*pm(Ir,Jr)
                    track(iygrd,itimep1,fid)=track(iygrd,itime,fid)+    &
     &                                  d_ij*sin(theta_ij)*pn(Ir,Jr)
                  END IF
                  track(izgrd,itimep1,fid)=REAL(KrMax,r8)
!
! -------- KINESIS-BASED SWIMMING BEHAVIOR (HUMSTON) --------
                ELSE IF ((swimtype(1,fid).ge.2).and.                    &
     &                   (swimtype(2,fid).ge.2)) THEN
!
! migration time for sardines
                  mg_time=REAL(INT(time(ng)/86400.0_r8/days_year))
                  mg_time=time(ng)/86400.0_r8-days_year*mg_time
!
                  h1=0.75_r8
                  h2=0.9_r8
                  tau=5.0_r8
! Temprature cue
                  fTemp=track(itemp+NFV(ng)-NT(ng),itime,fid)
                  fTemp_opt=0.5_r8*(FspTmin(isp,ng)+FspTmax(isp,ng))
                  fTemp_sig=2.0_r8
                  facT=EXP(-0.5_r8*((fTemp-fTemp_opt)/fTemp_sig)**2)
                  fTemp_alpha=1.0_r8+tau*(1.0_r8-facT)
! P-value cue
                  fPval=bioenergy(ifpval,fid)
                  fPval_opt=0.8_r8
                  fPval_sig=0.3_r8
                  fPval=MIN(fPval,fPval_opt)
                  facPv=EXP(-0.5_r8*((fPval-fPval_opt)/fPval_sig)**2)
                  fPval_alpha=1.0_r8+tau*(1.0_r8-facPv)
! Fish position (lon, lat)
                  zf_lon=track(iflon,itime,fid)
                  zf_lat=track(iflat,itime,fid)
!
! ----- Horizontal behavior (updated once daily) -----
!
                  IF (lifestage(fid).ge.if_juvenile) THEN

                    IF (swimtype(1,fid).eq.2) THEN
! Horizontal behavior based on temperature
                      cue1=facT
                      cue2=0.0_r8
                      cue1_a=fTemp_alpha
                      cue2_a=0.0_r8
                    ELSE IF (swimtype(1,fid).eq.3) THEN
! Horizontal behavior based on p-value
                      cue1=facPv
                      cue2=0.0_r8
                      cue1_a=fPval_alpha
                      cue2_a=0.0_r8
                    ELSE IF (swimtype(1,fid).eq.4) THEN
! Horinztal behavior based on temperature and p-value
                      cue1=facT
                      cue2=facPv
                      cue1_a=fTemp_alpha
                      cue2_a=fPval_alpha
                    END IF
!
! Displacement from previous time step
                    xdis=(track(ixgrd,itime,fid)-                       &
     &                    track(ixgrd,itimem1,fid))/pm(Ir,Jr)
                    ydis=(track(iygrd,itime,fid)-                       &
     &                    track(iygrd,itimem1,fid))/pn(Ir,Jr)
! Use xdis, ydis below for purely random motion
!                    xdis=0.0_r8
!                    ydis=0.0_r8
!
                    IF ((MOD(iic(ng)-1,24*steps_per_hour(ng)).eq.0).and.(iic(ng).ne.ntstart(ng))) THEN 
                    !IF (MOD(iic(ng)-1,24*steps_per_hour(ng)).eq.0) THEN ! RD bug
! Uswim is in cm/s, so multiply by 0.01 for m/s)
                      dist=0.01_r8*Uswim*dtfish
! Compute random displacements
                      epsd=(0.5_r8*dist**2)**0.5_r8
                      epsd_sigma=0.5_r8*dist
                      lflag=.TRUE.
                      DO WHILE (lflag)
                        CALL gasdev(snudg)
                        IF (ABS(snudg).le.3.0_r8) lflag=.FALSE.
                      END DO
                      ee=snudg*epsd_sigma+epsd
!                    ee=rwalk(fid)*epsd_sigma+epsd
!                    IF (r2walk(fid).lt.0.5_r8) THEN
!                      epsx=-ee
!                    ELSE
!                      epsx=ee
!                    END IF
!                    ee=rwalk(fid+Nfish(ng))*epsd_sigma+epsd
!                    IF (r2walk(fid).lt.0.5_r8) THEN
!                      epsy=-ee
!                    ELSE
!                      epsy=ee
!                    END IF
                      CALL ran1 (snudg)
# if defined WC12 || defined WC13 || defined CCS1
                      IF (isp.eq.if_anchovy) THEN
! Biased random mvt toward coast if offshore of 2000m isobath.
                        IF (h(Ir,Jr).le.2000.0_r8) THEN
                          IF (snudg.lt.0.5_r8) epsx=-ee
                          IF (snudg.ge.0.5_r8) epsx=ee
                        ELSE
                          IF (snudg.lt.0.3_r8) epsx=-ee
                          IF (snudg.ge.0.3_r8) epsx=ee
                        END IF
                      ELSE IF (isp.eq.if_sardine) THEN
                        IF (snudg.lt.0.5_r8) epsx=-ee
                        IF (snudg.ge.0.5_r8) epsx=ee
                      END IF
                      lflag=.TRUE.
                      DO WHILE (lflag)
                        CALL gasdev(snudg)
                        IF (ABS(snudg).le.3.0_r8) lflag=.FALSE.
                      END DO
                      ee=snudg*epsd_sigma+epsd
                      CALL ran1 (snudg)
                      IF (isp.eq.if_anchovy) THEN
! Biased random mvt toward coast if offshore of 2000m isobath.
                        IF (h(Ir,Jr).le.2000.0_r8) THEN
                          IF (snudg.lt.0.5_r8) epsy=-ee
                          IF (snudg.ge.0.5_r8) epsy=ee
                        ELSE
                          IF (zf_lat.le.39.5_r8) THEN
                            IF (snudg.lt.0.4_r8) epsy=-ee
                            IF (snudg.ge.0.4_r8) epsy=ee
                          ELSE IF (zf_lat.ge.40.5_r8) THEN
                            IF (snudg.lt.0.5_r8) epsy=-ee
                            IF (snudg.ge.0.5_r8) epsy=ee
                          ELSE
                            fac=0.4_r8+0.1_r8*(zf_lat-39.5_r8)
                            IF (snudg.lt.fac) epsy=-ee
                            IF (snudg.ge.fac) epsy=ee
                          END IF
                        END IF
                      ELSE IF (isp.eq.if_sardine) THEN
! seasonal migration for adult sardine: north in June, south in Sept.
                        IF (lifestage(fid).ge.if_subadult) THEN
                          IF ((mg_time.gt.150.0_r8).and.                &
     &                        (mg_time.lt.180.0_r8).and.                &
     &                        (zf_lat.lt.46.5_r8)) THEN
                            IF (snudg.lt.0.15_r8) epsy=-ee
                            IF (snudg.ge.0.15_r8) THEN
                              epsy=ee
                              IF ((zf_lon.ge.-124.5).and.                 &
     &                            (zf_lat.le.40.5_r8)) epsx=-0.7_r8*epsy
                            ENDIF
                          ELSE IF ((mg_time.gt.240.0_r8).and.           &
     &                             (mg_time.lt.270.0_r8).and.           &
     &                             (zf_lat.gt.31.5_r8)) THEN
                            IF (snudg.ge.0.15_r8) THEN
                              epsy=-ee
                              IF (zf_lat.le.40.5_r8) epsx=-0.7_r8*epsy
                            ENDIF
                            IF (snudg.lt.0.15_r8) epsy=ee
                          ELSE
                            IF (snudg.lt.0.5_r8) epsy=-ee
                            IF (snudg.ge.0.5_r8) epsy=ee
                          END IF
                        ELSE
                          IF (snudg.lt.0.5_r8) epsy=-ee
                          IF (snudg.ge.0.5_r8) epsy=ee
                        END IF
                      END IF
# endif
# if defined NWPACIFIC
! Shinichi custom code for nw pacific
                      IF (isp.eq.if_anchovy) THEN
                        IF (h(Ir,Jr).le.2000.0_r8) THEN
                          IF (snudg.lt.0.5_r8) epsx=-ee
                          IF (snudg.ge.0.5_r8) epsx=ee
                        ELSE
                          IF (snudg.lt.0.5_r8) epsx=-ee
                          IF (snudg.ge.0.5_r8) epsx=ee
                        END IF
                      ELSE IF (isp.eq.if_sardine) THEN
                        IF (snudg.lt.0.5_r8) epsx=-ee
                        IF (snudg.ge.0.5_r8) epsx=ee
                      END IF
                      lflag=.TRUE.
                      DO WHILE (lflag)
                        CALL gasdev(snudg)
                        IF (ABS(snudg).le.3.0_r8) lflag=.FALSE.
                      END DO
                      ee=snudg*epsd_sigma+epsd
                      CALL ran1 (snudg)
                      IF (isp.eq.if_anchovy) THEN
! Biased random mvt toward coast if offshore of 2000m isobath.
                        IF (h(Ir,Jr).le.2000.0_r8) THEN
                          IF (snudg.lt.0.5_r8) epsy=-ee
                          IF (snudg.ge.0.5_r8) epsy=ee
                        ELSE
                          IF (zf_lat.le.30.0_r8) THEN
                            IF (snudg.lt.0.4_r8) epsy=-ee
                            IF (snudg.ge.0.4_r8) epsy=ee
                          ELSE IF (zf_lat.ge.33.0_r8) THEN
                            IF (snudg.lt.0.5_r8) epsy=-ee
                            IF (snudg.ge.0.5_r8) epsy=ee
                          ELSE
                            fac=0.4_r8+0.1_r8*(zf_lat-30.0_r8)
                            IF (snudg.lt.fac) epsy=-ee
                            IF (snudg.ge.fac) epsy=ee
                          END IF
                        END IF
                      ELSE IF (isp.eq.if_sardine) THEN
! seasonal migration for adult sardine: north in June, south in Sept.
                        IF (lifestage(fid).ge.if_subadult) THEN
                          IF ((mg_time.gt.150.0_r8).and.                &
     &                        (mg_time.lt.180.0_r8).and.                &
     &                        (zf_lat.lt.46.5_r8)) THEN
                            IF (snudg.lt.0.15_r8) epsy=-ee
                            IF (snudg.ge.0.15_r8) THEN
                              epsy=ee
                              IF ((zf_lon.ge.140.0_r8).and.             &
     &                            (zf_lat.le.34.0_r8)) epsx=+0.7_r8*epsy
                            ENDIF
                          ELSE IF ((mg_time.gt.240.0_r8).and.           &
     &                             (mg_time.lt.270.0_r8).and.           &
     &                             (zf_lat.gt.33.0_r8)) THEN
                            IF (snudg.ge.0.15_r8) THEN
                              epsy=-ee
                              IF (zf_lat.le.46.0_r8) epsx=+0.7_r8*epsy
                            ENDIF
                            IF (snudg.lt.0.15_r8) epsy=ee
                          ELSE
                            IF (snudg.lt.0.5_r8) epsy=-ee
                            IF (snudg.ge.0.5_r8) epsy=ee
                          END IF
                        ELSE
                          IF (snudg.lt.0.5_r8) epsy=-ee
                          IF (snudg.ge.0.5_r8) epsy=ee
                        END IF
                      END IF
# endif
! Compute happiness-based  displacements
                      cue1_f=h1*cue1
                      cue1_g=1.0_r8-h2*cue1
                      cue2_f=h1*cue2
                      cue2_g=1.0_r8-h2*cue2
                      fac=cue1_a+cue2_a
                      f_xdis=xdis*(cue1_a*cue1_f+cue2_a*cue2_f)/fac
                      g_epsx=epsx*(cue1_a*cue1_g+cue2_a*cue2_g)/fac
                      f_ydis=ydis*(cue1_a*cue1_f+cue2_a*cue2_f)/fac
                      g_epsy=epsy*(cue1_a*cue1_g+cue2_a*cue2_g)/fac
                      xdis=f_xdis+g_epsx
                      ydis=f_ydis+g_epsy
                    END IF
! Swim towards best location
                    track(ixgrd,itimep1,fid)=track(ixgrd,itime,fid)+    &
     &                                                xdis*pm(Ir,Jr)
                    track(iygrd,itimep1,fid)=track(iygrd,itime,fid)+    &
     &                                                ydis*pn(Ir,Jr)
                  END IF
!
! ----- Vertical behavior (update at every time step) -----
!
                  IF (lifestage(fid).ge.if_juvenile) THEN
!
                    IF (swimtype(2,fid).eq.2) THEN
! Vertical behavior based on temperature
                      cue1=facT
                      cue2=0.0_r8
                      cue1_a=fTemp_alpha
                      cue2_a=0.0_r8
                    ELSE IF (swimtype(2,fid).eq.3) THEN
! Vertical behavior based on p-value
                      cue1=facPv
                      cue2=0.0_r8
                      cue1_a=fPval_alpha
                      cue2_a=0.0_r8
                    ELSE IF (swimtype(2,fid).eq.4) THEN
! Vertical behavior based on temperature and p-value
                      cue1=facT
                      cue2=facPv
                      cue1_a=fTemp_alpha
                      cue2_a=fPval_alpha
                    END IF
!
! Use commented lines below to check behavior based on depth
!                    Fvar=track(idpth,itime,fid)
!                    Fvar_opt=300.0_r8
!                    Fvar_sigma=100.0_r8
!
                    zdis=track(idpth,itime,fid)-track(idpth,itimem1,fid)
! Use zdis below for purely random motion
!                    zdis=0.0_r8
! Uswim is in cm/s, so multiply by 0.01 for m/s)
!                    dist=0.01_r8*Uswim*dtfish
! Use 10x cell thickness for random displacement
                    dist=10.0_r8*Hz(Ir,Jr,Kr)
! Compute random displacements
!                  epsd=(0.5_r8*dist**2)**0.5_r8
!                  epsd_sigma=0.5_r8*dist
!                  ee=rwalk(fid+2*Nfish(ng))*epsd_sigma+epsd
!                  IF (r2walk(fid).lt.0.5_r8) THEN
!                    epsz=-ee
!                  ELSE
!                    epsz=ee
!                  END IF
                    epsd=(0.5_r8*dist**2)**0.5_r8
                    epsd_sigma=0.5_r8*dist
                    lflag=.TRUE.
                    DO WHILE (lflag)
                      CALL gasdev(snudg)
                      IF (ABS(snudg).le.3.0_r8) lflag=.FALSE.
                    END DO
                    ee=snudg*epsd_sigma+epsd
                    CALL ran1 (snudg)
                    IF (snudg.lt.0.5_r8) epsz=-ee
                    IF (snudg.ge.0.5_r8) epsz=ee
! Compute happiness-based  displacements
                    cue1_f=h1*cue1
                    cue1_g=1.0_r8-h2*cue1
                    cue2_f=h1*cue2
                    cue2_g=1.0_r8-h2*cue2
                    fac=cue1_a+cue2_a
                    f_zdis=zdis*(cue1_a*cue1_f+cue2_a*cue2_f)/fac
                    g_epsz=epsz*(cue1_a*cue1_g+cue2_a*cue2_g)/fac
                    zdis=f_zdis+g_epsz
! Swim towards best location
                    zkdis=0.0_r8
                    KrMax=Kr
                    IF (zdis.gt.0.0_r8) THEN
                      DO kk=MIN(Kr,N(ng)),N(ng)
                        zkdis=zkdis+Hz(Ir,Jr,kk)
                        IF (zkdis.lt.zdis) KrMax=kk
                      END DO
                    END IF
                    IF (zdis.lt.0.0_r8) THEN
                      DO kk=MAX(Kr,2),2,-1
                        zkdis=zkdis+Hz(Ir,Jr,kk)
                        IF (zkdis.lt.ABS(zdis)) KrMax=kk
                      END DO
                    END IF
                    KrMax=MIN(N(ng),MAX(2,KrMax))
                    track(izgrd,itimep1,fid)=REAL(KrMax,r8)
! Impose 2D movement at surface
!                    track(izgrd,itimep1,fid)=REAL(N(ng),r8)
!
                  ELSE IF (lifestage(fid).eq.if_larva) THEN
                    scoreMax=0.0_r8
                    KrMax=N(ng)
                    DO kk=1,N(ng)
# ifdef MASKING
                      Bsz=MAX(0.0_r8,OCEAN(ng)%t(Ir,Jr,kk,nnew,iSzoo)*  &
     &                               rmask(Ir,Jr))
                      Blz=MAX(0.0_r8,OCEAN(ng)%t(Ir,Jr,kk,nnew,iLzoo)*  &
     &                               rmask(Ir,Jr))
                      Bpz=MAX(0.0_r8,OCEAN(ng)%t(Ir,Jr,kk,nnew,iPzoo)*  &
     &                               rmask(Ir,Jr))
# else
                      Bsz=MAX(0.0_r8,OCEAN(ng)%t(Ir,Jr,kk,nnew,iSzoo))
                      Blz=MAX(0.0_r8,OCEAN(ng)%t(Ir,Jr,kk,nnew,iLzoo))
                      Bpz=MAX(0.0_r8,OCEAN(ng)%t(Ir,Jr,kk,nnew,iPzoo))
# endif
                      Csz=Bsz*ZSpref(ils,isp,ng)/K_ZS(ils,isp,ng)
                      Clz=Blz*ZLpref(ils,isp,ng)/K_ZL(ils,isp,ng)
                      Cpz=Bpz*ZPpref(ils,isp,ng)/K_ZP(ils,isp,ng)
                      score=Csz+Clz+Cpz
                      IF (score.gt.scoreMax) THEN
                        scoreMax=score
                        KrMax=kk
                      END IF
                    END DO
                    track(izgrd,itimep1,fid)=REAL(KrMax,r8)
                  END IF
!
                END IF
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE fish_swim
#endif
      END MODULE fish_swim_mod
