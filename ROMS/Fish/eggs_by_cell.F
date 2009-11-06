#include "cppdefs.h"
      MODULE eggs_by_cell_mod
#if defined NEMURO_SAN && defined EGGS_BISECTION
!
!svn $Id: wetdry.F 975 2009-05-05 22:51:13Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!===================================================== Kate Hedstrom ===
!                                                                      !
!  This routine counts the number of fish per cell and builds a        !
!  linked list of them.                                                !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: eggs_by_cell

      CONTAINS
!
!***********************************************************************
      SUBROUTINE eggs_by_cell (ng, tile)

      USE mod_param
      USE mod_types
      USE mod_ocean
      USE mod_fish
#ifdef MASKING
      USE mod_grid
#endif
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, tile

# include "tile.h"

      CALL eggs_by_cell_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       OCEAN(ng) % egg_count,                     &
     &                       FISHES(ng) % bioenergy,                    &
     &                       FISHES(ng) % bounded,                      &
     &                       FISHES(ng) % alive,                        &
     &                       FISHES(ng) % species,                      &
     &                       FISHES(ng) % track)

      END SUBROUTINE eggs_by_cell

!***********************************************************************
      SUBROUTINE eggs_by_cell_tile (ng, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        rmask,                                    &
# endif
     &                        egg_count, bioenergy,                     &
     &                        bounded, alive, species, track)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      integer, intent(out) :: egg_count(LBi:,LBj:,:)
      real(r8), intent(inout) :: bioenergy(:,:)
      logical, intent(in) :: bounded(:)
      logical, intent(in) :: alive(:)
      integer, intent(in) :: species(:)
      real(r8), intent(in) :: track(:,0:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      integer, intent(out) :: egg_count(LBi:UBi,LBj:UBj,Nspecies(ng))
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfish(ng))
      logical, intent(in) :: bounded(Nfish(ng))
      logical, intent(in) :: alive(Nfish(ng))
      integer, intent(in) :: species(Nfish(ng))
      real(r8), intent(in) :: track(NFV(ng),0:NFT,Nfish(ng))
# endif
!
      integer :: i, j, l, isp
      integer :: Imin, Imax, Jmin, Jmax

      real(r8) :: cff

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
! Clean out to start. 
!-----------------------------------------------------------------------
!
      DO isp=1,Nspecies(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            egg_count(i,j,isp) = 0
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
! Loop over all fish and find the ones on my tile.
!-----------------------------------------------------------------------
!
      DO l=1,Nfish(ng)
        IF (bounded(l) .and. alive(l)) THEN
          i = int(track(ixgrd,nf(ng),l)+0.5)
          j = int(track(iygrd,nf(ng),l)+0.5)
          isp = idfish(species(l))
! bounded has been collected and no longer reflects whether the float is
! on this tile
          IF (i .ge. Imin .and. i .le. Imax .and.                       &
     &        j .ge. Jmin .and. j .le. Jmax) THEN
            egg_count(i,j,isp) = egg_count(i,j,isp) + bioenergy(ifeggs,l)
          END IF
        END IF
      END DO
!
      RETURN
      END SUBROUTINE eggs_by_cell_tile
#endif
      END MODULE eggs_by_cell_mod
