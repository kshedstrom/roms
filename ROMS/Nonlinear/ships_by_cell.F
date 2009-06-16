#include "cppdefs.h"
      MODULE ships_by_cell_mod
#if defined NEMURO_SAN && defined FLEET
!
!svn $Id: wetdry.F 975 2009-05-05 22:51:13Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!===================================================== Kate Hedstrom ===
!                                                                      !
!  This routine counts the number of ships per cell and builds a        !
!  linked list of them.                                                !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: ships_by_cell

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ships_by_cell (ng, tile)

      USE mod_param
      USE mod_types
      USE mod_ocean
      USE mod_ships
#ifdef MASKING
      USE mod_grid
#endif
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, tile

# include "tile.h"

      CALL ships_by_cell_tile (ng, tile,                                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       OCEAN(ng) % ships_count,                   &
     &                       OCEAN(ng) % ships_list,                    &
     &                       SHIPS(ng) % shipnodes,                     &
     &                       SHIPS(ng) % bounded,                       &
     &                       SHIPS(ng) % track)

      END SUBROUTINE ships_by_cell

!***********************************************************************
      SUBROUTINE ships_by_cell_tile (ng, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        rmask,                                    &
# endif
     &                        ships_count, ships_list, shipnodes,         &
     &                        bounded, track)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_types
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
      integer, intent(out) :: ships_count(LBi:,LBj:)
      type(fishnode), intent(out) :: ships_list(LBi:,LBj:)
      type(fishnode), target, intent(inout) :: shipnodes(:)
      logical, intent(in) :: bounded(:)
      real(r8), intent(in) :: track(:,0:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      integer, intent(out) :: ships_count(LBi:UBi,LBj:UBj)
      type(fishnode), intent(out) :: ships_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(inout) :: shipnodes(Nfloats(ng))
      logical, intent(in) :: bounded(Nfloats(ng))
      real(r8), intent(in) :: track(NFV(ng),0:NFT,Nfloats(ng))
# endif
!
      integer :: i, j, l
      integer :: Imin, Imax, Jmin, Jmax

      real(r8) :: cff

# include "set_bounds.h"
# ifdef DISTRIBUTED
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
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ships_count(i,j) = 0
          ships_list(i,j) % next => null()
        END DO
      END DO
!
!-----------------------------------------------------------------------
! Loop over all ships.
!-----------------------------------------------------------------------
!
      DO l=1,Nships(ng)
        IF (bounded(l)) THEN
          i = int(track(ixgrd,nf(ng),l)+0.5)
          j = int(track(iygrd,nf(ng),l)+0.5)
! bounded has been collected and no longer reflects whether the float is
! on this tile
          IF (i .ge. Imin .and. i .le. Imax .and.                       &
     &        j .ge. Jmin .and. j .le. Jmax) THEN
            ships_count(i,j) = ships_count(i,j) + 1
            shipnodes(l) % next => ships_list(i,j) % next
            ships_list(i,j) % next => shipnodes(l)
          END IF
        END IF
      END DO
!
      RETURN
      END SUBROUTINE ships_by_cell_tile
#endif
      END MODULE ships_by_cell_mod
