#include "cppdefs.h"
      MODULE fish_by_cell_mod
#ifdef NEMURO_SAN
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
      PUBLIC :: fish_by_cell

      CONTAINS
!
!***********************************************************************
      SUBROUTINE fish_by_cell (ng, tile)

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

      CALL fish_by_cell_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       OCEAN(ng) % fish_count,                    &
     &                       OCEAN(ng) % fish_list,                     &
     &                       FISHES(ng) % fishnodes,                    &
     &                       OCEAN(ng) % fish0,                         &
     &                       FISHES(ng) % bounded,                      &
     &                       FISHES(ng) % alive,                        &
     &                       FISHES(ng) % track)

      END SUBROUTINE fish_by_cell

!***********************************************************************
      SUBROUTINE fish_by_cell_tile (ng, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        rmask,                                    &
# endif
     &                        fish_count, fish_list, fishnodes,         &
     &                        fish0, bounded, alive, track)
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
      type(fishnode), target, intent(in) :: fish0
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      integer, intent(out) :: fish_count(LBi:,LBj:)
      type(fishnode), intent(out) :: fish_list(LBi:,LBj:)
      type(fishnode), target, intent(inout) :: fishnodes(:)
      logical, intent(in) :: bounded(:)
      logical, intent(in) :: alive(:)
      real(r8), intent(in) :: track(:,0:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      integer, intent(out) :: fish_count(LBi:UBi,LBj:UBj)
      type(fishnode), intent(out) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(inout) :: fishnodes(Nfloats(ng))
      logical, intent(in) :: bounded(Nfloats(ng))
      logical, intent(in) :: alive(Nfloats(ng))
      real(r8), intent(in) :: track(NFV(ng),0:NFT,Nfloats(ng))
# endif
!
      integer :: i, j, l

      real(r8) :: cff
      type(fishnode), pointer :: newfish

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
! Clean out to start. 
!-----------------------------------------------------------------------
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          fish_count(i,j) = 0
          fish_list(i,j) % next => fish0
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
          fish_count(i,j) = fish_count(i,j) + 1
          fishnodes(l) % next => fish_list(i,j) % next
          fish_list(i,j) % next => fishnodes(l)
        END IF
      END DO
!
      RETURN
      END SUBROUTINE fish_by_cell_tile
#endif
      END MODULE fish_by_cell_mod
