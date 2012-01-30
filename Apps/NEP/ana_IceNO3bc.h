      SUBROUTINE ana_IceNO3bc (ng, tile)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets free-surface open boundary conditions using       !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL ana_IceNO3bc_tile (ng, tile,                                   &
     &                     LBi, UBi, LBj, UBj)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(45)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_IceNO3bc
!
!***********************************************************************
      SUBROUTINE ana_IceNO3bc_tile (ng, tile,                             &
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: cff, fac, omega, phase, val

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Free-surface open boundary conditions.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%IceNO3_east(j)=0.0_r8
        END DO
      END IF
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrR,JendR
          BOUNDARY(ng)%IceNO3_west(j)=0.0_r8
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%IceNO3_south(i)=0.0_r8
        END DO
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrR,IendR
          BOUNDARY(ng)%IceNO3_north(i)=0.0_r8
        END DO
      END IF
      RETURN
      END SUBROUTINE ana_IceNO3bc_tile