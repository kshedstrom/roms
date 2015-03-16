#include "cppdefs.h"
      MODULE new_year_fleet_mod
#if defined NEMURO_SAN && defined FISHING_FLEET
!
!svn $Id$
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  At the end of each year, we kill off any remaining fish of the
!  the oldest age class to make room for the next age class. We also
!  move the indices for creating new fish.
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: new_year_fleet

      CONTAINS
!
!***********************************************************************
      SUBROUTINE new_year_fleet(ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_fish
      USE mod_fleet
      USE mod_fishing
      USE mod_ocean
      USE mod_grid
      USE mod_scalars
      USE mod_stepping
      USE mod_parallel
      USE mod_biology
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer :: ib, ifid
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif

      CALL gather_array(ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                 GRID(ng)%xr, GRID(ng)%yr, OCEAN(ng)%fish_count,  &
     &                 OCEAN(ng)%distport, OCEAN(ng)%nearport,          &
     &                 OCEAN(ng)%cpue, OCEAN(ng)%avgcpue,               &
                       OCEAN(ng)%percpue, OCEAN(ng)%new_climcpue,       &
                       OCEAN(ng)%climcpue)
!
      IF (Master) THEN
!
        BOATS(ng)%AnnCatch=0.0_r8
        BOATS(ng)%TotFishBio=0.0_r8
        DO ib=1,Nboats(ng)
          BOATS(ng)%cumcatch(ib)=0.0_r8
        END DO
        DO ifid=1,Nfish(ng)
          IF ((FISHES(ng)%species(ifid).eq.if_sardine).and.             &
     &        (FISHES(ng)%alive(ifid)).and.                             &
     &        (FISHES(ng)%lifestage(ifid).ge.if_subadult)) THEN
            BOATS(ng)%TotFishBio=BOATS(ng)%TotFishBio+                  &
     &                      0.001_r8*FISHES(ng)%bioenergy(ifwwt,ifid)*  &
     &                              FISHES(ng)%bioenergy(ifworth,ifid)
          END IF
        END DO
        CALL updatecpue_clim(ng)
!
      END IF
!
      CALL scatter_array(ng, iNLM, LBi, UBi, LBj, UBj,                  &
     &                   OCEAN(ng)%distport, OCEAN(ng)%nearport,        &
     &                   OCEAN(ng)%cpue, OCEAN(ng)%avgcpue,             &
     &                   OCEAN(ng)%percpue, OCEAN(ng)%new_climcpue,     &
     &                   OCEAN(ng)%climcpue)
!
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE new_year_fleet
!
#endif
      END MODULE new_year_fleet_mod
