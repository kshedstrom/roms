#include "cppdefs.h"
      MODULE ships_by_cell_mod
#ifdef NEMURO_SAN
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
     &                       OCEAN(ng) % fish0,                         &
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
     &                        fish0, bounded, track)
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
      type(fishnode), pointer, intent(in) :: fish0
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

      real(r8) :: cff

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
! Clean out to start. 
!-----------------------------------------------------------------------
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          ships_count(i,j) = 0
          ships_list(i,j) % next => fish0
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
          ships_count(i,j) = ships_count(i,j) + 1
          shipnodes(l) % next => ships_list(i,j) % next
          ships_list(i,j) % next => shipnodes(l)
        END IF
      END DO
!
      RETURN
      END SUBROUTINE ships_by_cell_tile
#endif
      END MODULE ships_by_cell_mod
