#include "cppdefs.h"
      MODULE fish_elh_mod
#if defined NONLINEAR && defined NEMURO_SAN
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group         Mark Hadfield   !
!    Licensed under a MIT/X style license             John M. Klinck   !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes early life history of fish.                   !
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: fish_elh

      CONTAINS
!
!***********************************************************************
      SUBROUTINE count_young (ng, tile, LBi, UBi, LBj, UBj,             &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        species, lifestage,                       &
     &                        juv_cell, fish_count, fish_list,          &
     &                        bioenergy)
!***********************************************************************
!
      USE mod_param
      USE mod_fish
      USE mod_biology, only: idfish, max_species, DDscale_j
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef ASSUMED_SHAPE
      integer, intent(in) :: species(:)
      integer, intent(in) :: lifestage(:)
      real(r8), intent(out) :: juv_cell(IminS:,JminS:,:)
      integer, intent(in) :: fish_count(LBi:,LBj:)
      type(fishnode), intent(in) :: fish_list(LBi:,LBj:)
      real(r8), intent(in) :: bioenergy(:,:)
# else
      integer, intent(in) :: species(Nfish(ng))
      integer, intent(in) :: lifestage(Nfish(ng))
      real(r8), intent(out) :: juv_cell(IminS:ImaxS,JminS:JmaxS,        &
     &                         max_species)
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)
      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bioenergy(NFishV(ng),Nfish(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i1, i2, j1, j2, i, j, ii, jj, k
      integer :: ifish, ifid, isp, ifsp, ils
      type(fishnode), pointer :: thisfish

# include "set_bounds.h"

      juv_cell = 0
      DO i=Istr,Iend
        DO j=Jstr,Jend
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              isp = idfish(species(ifid))
              ils=lifestage(ifid)
              IF (ils.le.if_juvenile) THEN
                juv_cell(i,j,isp) = juv_cell(i,j,isp) +                 &
     &               bioenergy(ifworth,ifid)
              END IF
              thisfish => thisfish % next
            END DO
          END IF
          DO isp=1,max_species
            juv_cell(i,j,isp) = MIN(juv_cell(i,j,isp)*DDscale_J(isp,ng),&
     &                    100._r8)
          END DO
        END DO
      END DO

      END SUBROUTINE count_young
!
!
!***********************************************************************
      SUBROUTINE fish_elh (ng, tile, LBi, UBi, LBj, UBj, LBk, UBk,      &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     itime, nnew,                                 &
     &                     fishthread, bounded, track,                  &
     &                     bioenergy, alive, species,                   &
     &                     lifestage, deathby, fmortN,                  &
     &                     egg_dur, egg_num,                            &
     &                     ysac_dur, ysac_num,                          &
     &                     larv_dur, larv_num,                          &
     &                     juv_dur, juv_num, suba_num,                  &
     &                     fish_count, fish_list, fishnodes)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
      USE mod_types
      USE mod_grid
      USE mod_biology
      USE mod_fish
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, itime, nnew
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef ASSUMED_SHAPE
      integer, intent(in) :: fish_count(LBi:,LBj:)
      integer, intent(inout) :: species(:)
      integer, intent(inout) :: lifestage(:)
      integer, intent(inout) :: deathby(:)

      type(fishnode), intent(in) :: fish_list(LBi:,LBj:)
      type(fishnode), target, intent(in) :: fishnodes(:)

      logical, intent(in) :: fishthread(:)
      logical, intent(inout) :: bounded(:)
      logical, intent(inout) :: alive(:)

      real(r8), intent(inout) :: track(:,0:,:)
      real(r8), intent(inout) :: bioenergy(:,:)
      real(r8), intent(inout) :: fmortN(:)

      real(r8), intent(inout) :: egg_dur(:)
      real(r8), intent(inout) :: egg_num(:)
      real(r8), intent(inout) :: ysac_dur(:)
      real(r8), intent(inout) :: ysac_num(:)
      real(r8), intent(inout) :: larv_dur(:)
      real(r8), intent(inout) :: larv_num(:)
      real(r8), intent(inout) :: juv_dur(:)
      real(r8), intent(inout) :: juv_num(:)
      real(r8), intent(inout) :: suba_num(:)
# else
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)
      integer, intent(inout) :: species(Nfish(ng))
      integer, intent(inout) :: lifestage(Nfish(ng))
      integer, intent(inout) :: deathby(Nfish(ng))

      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Nfish(ng))

      logical, intent(in) :: fishthread(Nfish(ng))
      logical, intent(inout) :: bounded(Nfish(ng))
      logical, intent(inout) :: alive(Nfish(ng))

      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfish(ng))
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfish(ng))
      real(r8), intent(inout) :: fmortN(Nfish(ng))

      real(r8), intent(inout) :: egg_dur(Nfish(ng))
      real(r8), intent(inout) :: egg_num(Nfish(ng))
      real(r8), intent(inout) :: ysac_dur(Nfish(ng))
      real(r8), intent(inout) :: ysac_num(Nfish(ng))
      real(r8), intent(inout) :: larv_dur(Nfish(ng))
      real(r8), intent(inout) :: larv_num(Nfish(ng))
      real(r8), intent(inout) :: juv_dur(Nfish(ng))
      real(r8), intent(inout) :: juv_num(Nfish(ng))
      real(r8), intent(inout) :: suba_num(Nfish(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i1, i2, j1, j2, i, j, ii, jj, k
      integer :: ifish, ifid, isp, ifsp, ils

      type(fishnode), pointer :: thisfish
      real(r8) :: juv_cell(IminS:ImaxS,JminS:JmaxS,max_species)

      real(r8) :: Fweight, Fworth, Flength, Ftemp
      real(r8) :: dtdays, dtsecs, dy_time, t_ramp, tday
      real(r8) :: fac, a_t, d_t, rx

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute growth and/or mortality terms for early life stages
!-----------------------------------------------------------------------
!
!!      dtdays=dt(ng)*sec2day
!      dtsecs=4.0_r8*dt(ng)
!      dtdays=4.0_r8*dt(ng)*sec2day
! RD BUG above, we call fish every hour dtsecs = 3600, dtdays = 3600 * sec2day
! if we decide to change the frequency we call fish, we'll have to change this :
      dtsecs=3600.0_r8
      dtdays=3600.0_r8*sec2day

! dy_time: time in days, modulo days_year
! Note: only works if time origin is at the start of a year
      dy_time=REAL(INT(time(ng)/86400.0_r8/days_year))
      dy_time=time(ng)/86400.0_r8-days_year*dy_time
      CALL count_young (ng, tile, LBi, UBi, LBj, UBj,                   &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  species, lifestage,                             &
     &                  juv_cell, fish_count, fish_list, bioenergy)
!
!      DO i=LBi+2,UBi-2
!        DO j=LBj+2,UBj-2
      DO i=Istr,Iend
        DO j=Jstr,Jend
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              isp = idfish(species(ifid))
              ils=lifestage(ifid)
              IF (ils.le.if_juvenile) THEN
                Fweight=bioenergy(ifwwt,ifid)
                Flength=bioenergy(iflngth,ifid)
                Fworth=bioenergy(ifworth,ifid)
                Ftemp=track(itemp+NFV(ng)-NT(ng),itime,ifid)
! Juvenile: mortality only (growth is done in fish_growth.F)
                IF (ils.eq.if_juvenile) THEN
!  rx is density dependent term
                  rx = (DDmort1_J(isp,ng) +                            &
     &                  juv_cell(i,j,isp)*DDmort2_J(isp,ng)) /         &
     &                 (DDmort3_J(isp,ng) +                            &
     &                  juv_cell(i,j,isp)*DDmort4_J(isp,ng))
                  Fworth=Fworth*EXP(-Nmort(ils,isp,ng)*rx*dtdays)
                  fmortN(ifid)=Nmort(ils,isp,ng)
                  juv_dur(ifid)=juv_dur(ifid)+dtdays
! JF NOTE: For now, J->A is at 1st birthday (1 Jan)
!                  IF (Flength.ge.LengthJA(isp,ng)) THEN
!                    lifestage(ifid)=if_subadult
!                    suba_num(ifid)=bioenergy(ifworth,ifid)
!                  ELSE
!                    juv_dur(ifid)=juv_dur(ifid)+dtdays
!                  END IF
! Larva: mortality only (growth is done in fish_growth.F)
                ELSE IF (ils.eq.if_larva) THEN
                  Fworth=Fworth*EXP(-Nmort(ils,isp,ng)*dtdays)
                  fmortN(ifid)=Nmort(ils,isp,ng)
!                  IF ((Fweight.ge.WeightLJ(isp,ng)).or.                 &
!     &                (Flength.ge.LengthLJ(isp,ng))) THEN
                  IF (Flength.ge.LengthLJ(isp,ng)) THEN
                    lifestage(ifid)=if_juvenile
                    juv_dur(ifid)=0.0_r8
                    juv_num(ifid)=bioenergy(ifworth,ifid)
                  ELSE
                    larv_dur(ifid)=larv_dur(ifid)+dtdays
                  END IF
! Yolk Sac
                ELSE IF (ils.eq.if_yolksac) THEN
                  IF (isp.eq.if_anchovy) THEN
                    a_t=0.00115_r8*EXP(5.5338_r8*(1.0_r8-               &
     &                    EXP(-0.0546_r8*Ftemp)))
                    d_t=LOG(5.1493_r8/(5.1493_r8-                       &
     &                    LOG(4.2_r8/2.86_r8)))/a_t
                  ELSE IF (isp.eq.if_sardine) THEN
                    a_t=0.00102_r8*EXP(6.1876_r8*(1.0_r8-               &
     &                    EXP(-0.049_r8*Ftemp)))
                    d_t=LOG(5.1918_r8/(5.1918_r8-                       &
     &                    LOG(5.97_r8/3.74_r8)))/a_t
                  END IF
                  Fweight=Fweight+dtdays/d_t
                  Fworth=Fworth*EXP(-Nmort(ils,isp,ng)*dtdays)
                  fmortN(ifid)=Nmort(ils,isp,ng)
                  IF (Fweight.ge.1.0_r8) THEN
                    lifestage(ifid)=if_larva
                    Fweight=Wffeed(isp,ng)
                    Flength=Lffeed(isp,ng)
                    larv_dur(ifid)=0.0_r8
                    larv_num(ifid)=bioenergy(ifworth,ifid)
                  ELSE
                    ysac_dur(ifid)=ysac_dur(ifid)+dtdays
                  END IF
! Egg
                ELSE IF (ils.eq.if_egg) THEN
                  IF (isp.eq.if_anchovy) THEN
                    a_t=0.0412_r8*EXP(5.5338_r8*(1.0_r8-                &
     &                    EXP(-0.0546_r8*Ftemp)))
                    d_t=LOG(5.1493_r8/(5.1493_r8-                       &
     &                    LOG(2.86_r8/0.025_r8)))/a_t
                  ELSE IF (isp.eq.if_sardine) THEN
                    a_t=0.0317_r8*EXP(6.1876_r8*(1.0_r8-                &
     &                    EXP(-0.049_r8*Ftemp)))
                    d_t=LOG(5.1918_r8/(5.1918_r8-                       &
     &                    LOG(3.74_r8/0.0341_r8)))/a_t
                  END IF
                  Fweight=Fweight+dtdays/d_t
                  Fworth=Fworth*EXP(-Nmort(ils,isp,ng)*dtdays)
                  fmortN(ifid)=Nmort(ils,isp,ng)
                  IF (Fweight.ge.1.0_r8) THEN
                    lifestage(ifid)=if_yolksac
                    Fweight=0.0_r8
                    ysac_dur(ifid)=0.0_r8
                    ysac_num(ifid)=bioenergy(ifworth,ifid)
                  ELSE
                    egg_dur(ifid)=egg_dur(ifid)+dtdays
                  END IF
                END IF
                bioenergy(ifwwt,ifid)=Fweight
                bioenergy(iflngth,ifid)=Flength
                bioenergy(ifworth,ifid)=Fworth
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE fish_elh
#endif
      END MODULE fish_elh_mod
