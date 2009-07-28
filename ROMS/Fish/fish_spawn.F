#include "cppdefs.h"
      MODULE fish_spawn_mod
#ifdef NEMURO_SAN
!
!svn $Id: interp_floats.F 294 2009-01-09 21:37:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group         Mark Hadfield   !
!    Licensed under a MIT/X style license             John M. Klinck   !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes numbers and locations of eggs spawned today.  !
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: fish_spawn

      CONTAINS
!
!***********************************************************************
      SUBROUTINE fish_spawn (ng, tile)
!***********************************************************************

      USE mod_param
      USE mod_types
      USE mod_fish
      USE mod_ocean
      USE mod_grid
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
      CALL fish_spawn_tile (ng, tile, nnew(ng), LBi, UBi, LBj, UBj,     &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        GRID(ng) % spawn_dist,                    &
     &                        OCEAN(ng) % fish_count,                   &
     &                        OCEAN(ng) % fish_list,                    &
     &                        FISHES(ng) % fishnodes,                   &
     &                        FISHES(ng) % bounded,                     &
     &                        FISHES(ng) % track,                       &
     &                        FISHES(ng) % bioenergy,                   &
     &                        FISHES(ng) % alive,                       &
     &                        FISHES(ng) % species,                     &
     &                        FISHES(ng) % sp_start,                    &
     &                        FISHES(ng) % sp_end,                      &
     &                        FISHES(ng) % num_free,                    &
     &                        FISHES(ng) % num_super)
!
      END SUBROUTINE fish_spawn
!
!***********************************************************************
      SUBROUTINE fish_spawn_tile (ng, tile, nnew, LBi, UBi, LBj, UBj,   &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              spawn_dist,                         &
     &                              fish_count, fish_list, fishnodes,   &
     &                              bounded, track, bioenergy, alive,   &
     &                              species, sp_start, sp_end,          &
     &                              num_free, num_super)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_types
      USE mod_biology
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nnew
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: spawn_dist(LBi:,LBj:)
      integer, intent(in) :: fish_count(LBi:,LBj:)
      type(fishnode), intent(in) :: fish_list(LBi:,LBj:)
      type(fishnode), target, intent(in) :: fishnodes(:)
      logical, intent(in) :: bounded(:)
      real(r8), intent(in) :: track(:,0:,:)
      real(r8), intent(inout) :: bioenergy(:,:)
      logical, intent(in) :: alive(:)
      integer, intent(in) :: species(:)
      real(r8), intent(in) :: sp_start(:)
      real(r8), intent(in) :: sp_end(:)
      integer, intent(in) :: num_free(:)
      integer, intent(out) :: num_super(:)
# else
      real(r8), intent(in) :: spawn_dist(LBi:UBi,LBj:UBj)
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)
      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Nfish(ng))
      logical, intent(in) :: bounded(Nfish(ng))
      real(r8), intent(in) :: track(NFV(ng),0:NFT,Nfish(ng))
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfish(ng))
      logical, intent(in) :: alive(Nfish(ng))
      integer, intent(in) :: species(Nfish(ng))
      real(r8), intent(in) :: sp_start(Nspecies(ng))
      real(r8), intent(in) :: sp_end(Nspecies(ng))
      integer, intent(in) :: num_free(Nspecies(ng))
      integer, intent(out) :: num_super(Nspecies(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, isp, ifish, ifid
      integer :: Imin, Imax, Jmin, Jmax

      type(fishnode), pointer :: thisfish

      real(r8) :: Ftemp, temp_opt, temp_sigma
      real(r8) :: sp_time

# include "set_bounds.h"
# ifdef DISTRIBUTE
      Imin = LBi
      Imax = UBi
      Jmin = LBj
      Jmax = UBj
# else
      Imin = Istr
      Imax = Iend
      Jmin = Jstr
      Jmax = Jend
# endif
!
!-----------------------------------------------------------------------
!  Compute feeding terms on zooplankton to update growth
!-----------------------------------------------------------------------
!
      temp_opt=15.0_r8
      temp_sigma=3.0_r8
!
! sp_time: time in days, modulo 360 days
      sp_time=REAL(INT(time(ng)/86400.0_r8/360.0_r8))
      sp_time=time(ng)/86400.0_r8-360.0_r8*sp_time
!
! find out how many free superindividuals per species for today
      DO isp=1, Nspecies(ng)
        num_super(isp) = num_free(isp) / (sp_end(isp) - sp_time + 1)
      END DO
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              IF (bounded(ifid).and.alive(ifid)) THEN
                isp = idfish(species(ifid))
! If inside spawning window compute spawning attributes
                IF ((sp_time.ge.sp_start(isp)).and.                     &
     &              (sp_time.le.sp_end(isp))) THEN
                  Ftemp=track(itemp+NFV(ng)-NT(ng),nnew,ifid)
                  IF ((alive(ifid)).and.                                &
     &              (bioenergy(ifage,ifid).ge.Fage0(ng)).and.           &
     &              (INT(bioenergy(ifbatch,ifid)).lt.Nbatch(ng)).and.   &
     &              (ABS(Ftemp-temp_opt).le.temp_sigma).and.            &
     &              ((sp_time-bioenergy(iftspwn,ifid)).ge.7.0_r8)) THEN
                    bioenergy(ifeggs,ifid)=20000.0_r8*                  &
     &                                     bioenergy(ifworth,ifid)
                    bioenergy(ifspwnloc,ifid)=spawn_dist(i,j)
                    bioenergy(ifbatch,ifid)=bioenergy(ifbatch,ifid)+    &
     &                                      1.0_r8
                    bioenergy(iftspwn,ifid)=sp_time
                  ELSE
                    bioenergy(ifeggs,ifid)=0.0_r8
                    bioenergy(ifspwnloc,ifid)=0.0_r8
                  END IF
                END IF
! If after end of spawning window, zero out spawning attributes
                IF (sp_time.gt.sp_end(isp)) THEN
                  bioenergy(ifeggs,ifid)=0.0_r8
                  bioenergy(ifspwnloc,ifid)=0.0_r8
                END IF
!
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE fish_spawn_tile
#endif
      END MODULE fish_spawn_mod
