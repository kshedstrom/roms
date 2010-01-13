#include "cppdefs.h"
      MODULE pred_by_cell_mod
#if defined NEMURO_SAN && defined PREDATORS
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!===================================================== Kate Hedstrom ===
!                                                                      !
!  This routine counts the number of pred per cell and builds a        !
!  linked list of them.                                                !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: pred_by_cell

      CONTAINS
!
!***********************************************************************
      SUBROUTINE pred_by_cell (ng, tile)

      USE mod_param
      USE mod_types
      USE mod_ocean
      USE mod_predators
#ifdef MASKING
      USE mod_grid
#endif
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, tile

# include "tile.h"

      CALL pred_by_cell_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       OCEAN(ng) % pred_count,                    &
     &                       OCEAN(ng) % pred_list,                     &
     &                       PREDATORS(ng) % shipnodes,                 &
     &                       PREDATORS(ng) % bounded,                   &
     &                       PREDATORS(ng) % track)

      END SUBROUTINE pred_by_cell

!***********************************************************************
      SUBROUTINE pred_by_cell_tile (ng, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        rmask,                                    &
# endif
     &                        pred_count, pred_list, shipnodes,         &
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
      integer, intent(out) :: pred_count(LBi:,LBj:)
      type(fishnode), intent(out) :: pred_list(LBi:,LBj:)
      type(fishnode), target, intent(inout) :: shipnodes(:)
      logical, intent(in) :: bounded(:)
      real(r8), intent(in) :: track(:,0:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      integer, intent(out) :: pred_count(LBi:UBi,LBj:UBj)
      type(fishnode), intent(out) :: pred_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(inout) :: shipnodes(Nfloats(ng))
      logical, intent(in) :: bounded(Nfloats(ng))
      real(r8), intent(in) :: track(NFV(ng),0:NFT,Nfloats(ng))
# endif
!
      integer :: i, j, l
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
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          pred_count(i,j) = 0
          pred_list(i,j) % next => null()
        END DO
      END DO
!
!-----------------------------------------------------------------------
! Loop over all pred.
!-----------------------------------------------------------------------
!
      DO l=1,Npred(ng)
        IF (bounded(l)) THEN
          i = int(track(ixgrd,nf(ng),l)+0.5)
          j = int(track(iygrd,nf(ng),l)+0.5)
! bounded has been collected and no longer reflects whether the float is
! on this tile
          IF (i .ge. Imin .and. i .le. Imax .and.                       &
     &        j .ge. Jmin .and. j .le. Jmax) THEN
            pred_count(i,j) = pred_count(i,j) + 1
            prednodes(l) % next => pred_list(i,j) % next
            pred_list(i,j) % next => prednodes(l)
          END IF
        END IF
      END DO
!
      RETURN
      END SUBROUTINE pred_by_cell_tile
#endif
      END MODULE pred_by_cell_mod
