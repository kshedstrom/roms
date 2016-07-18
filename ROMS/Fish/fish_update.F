#include "cppdefs.h"
      MODULE fish_update_mod
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
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: fish_update

      CONTAINS
!
!***********************************************************************
      SUBROUTINE fish_update (ng, tile, LBi, UBi, LBj, UBj, LBk, UBk,   &
     &                            itime, nnew,                          &
     &                            fishthread, bounded, track,             &
     &                            bioenergy, alive, species,            &
     &                            lifestage, deathby,                   &
     &                            egg_dur, egg_num,                     &
     &                            ysac_dur, ysac_num,                   &
     &                            larv_dur, larv_num,                   &
     &                            juv_dur, juv_num, suba_num,           &
     &                            fish_count, fish_list, fishnodes)
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
      real(r8), intent(inout) :: egg_dur(Nfish(ng))
      real(r8), intent(inout) :: egg_num(Nfish(ng))
      real(r8), intent(inout) :: ysac_dur(Nfish(ng))
      real(r8), intent(inout) :: ysac_num(Nfish(ng))
      real(r8), intent(inout) :: larv_dur(Nfish(ng))
      real(r8), intent(inout) :: larv_num(Nfish(ng))
      real(r8), intent(inout) :: juv_dur(Nfish(ng))
      real(r8), intent(inout) :: juv_num(Nfish(ng))
      real(r8), intent(inout) :: suba_num(Nfish(ng))
!
!  Local variable declarations.
!
      integer :: i1, i2, j1, j2, i, j, ii, jj, k
      integer :: ifish, ifid, isp, ifsp, ils

      type(fishnode), pointer :: thisfish

      real(r8) :: Fweight, Fworth, Flength, Ftemp
      real(r8) :: dtdays, dtsecs, dy_time, t_ramp, tday

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Update fish attributes
!-----------------------------------------------------------------------
!
!      dtdays=dt(ng)*sec2day
!      dtsecs=4.0_r8*dt(ng)
!      dtdays=4.0_r8*dt(ng)*sec2day
      dtsecs=3600.0_r8
      dtdays=3600.0_r8*sec2day
! dy_time: time in days, modulo days_year
! Note: only works if time origin is at the start of a year
      dy_time=REAL(INT(time(ng)/86400.0_r8/days_year))
      dy_time=time(ng)/86400.0_r8-days_year*dy_time
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
              IF (fishthread(ifid).and.bounded(ifid)) THEN
! From Kenny's code:
                IF (bioenergy(ifworth,ifid).lt.1.0_r8) THEN
                  alive(ifid)=.FALSE.
                  IF ((lifestage(ifid).le.if_juvenile).and.             &
     &                (deathby(ifid).eq.0)) THEN
                    deathby(ifid)=2
                  END IF
                END IF
                IF (alive(ifid)) THEN
! Update fish age
                  bioenergy(ifage,ifid)=bioenergy(ifage,ifid)+          &
     &                                        dtdays/days_year
                ELSE
                  IF (lifestage(ifid).eq.if_juvenile) THEN
                    juv_num(ifid)=0.0_r8
                    juv_dur(ifid)=0.0_r8
                  ELSE IF (lifestage(ifid).eq.if_larva) THEN
                    larv_num(ifid)=0.0_r8
                    larv_dur(ifid)=0.0_r8
                  ELSE IF (lifestage(ifid).eq.if_yolksac) THEN
                    ysac_num(ifid)=0.0_r8
                    ysac_dur(ifid)=0.0_r8
                  ELSE IF (lifestage(ifid).eq.if_egg) THEN
                    egg_num(ifid)=0.0_r8
                    egg_dur(ifid)=0.0_r8
                  END IF
!                  bioenergy(ifwwt,ifid)=0.0_r8
!                  bioenergy(ifworth,ifid)=0.0_r8
!                  bioenergy(ifage,ifid)=0.0_r8
!                  bioenergy(ifbday,ifid)=0.0_r8
!                  lifestage(ifid)=0
                  bounded(ifid) = .FALSE.
                END IF
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE fish_update
#endif
      END MODULE fish_update_mod
