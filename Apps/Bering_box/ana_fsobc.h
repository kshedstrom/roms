      SUBROUTINE ana_fsobc (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
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
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_fsobc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 6)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_fsobc
!
!***********************************************************************
      SUBROUTINE ana_fsobc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
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
# ifdef SOUTH_FSOBC
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrT,IendT
          BOUNDARY(ng)%zeta_south(i)=0.0_r8
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE ana_fsobc_tile
