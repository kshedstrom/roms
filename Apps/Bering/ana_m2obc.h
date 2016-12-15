      SUBROUTINE ana_m2obc (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets 2D momentum open boundary conditions using        !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_m2obc_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     knew(ng),                                    &
     &                     GRID(ng) % angler,                           &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % on_u,                             &
     &                     GRID(ng) % om_v,                             &
#ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#endif
     &                     OCEAN(ng) % zeta)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(12)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_m2obc
!
!***********************************************************************
      SUBROUTINE ana_m2obc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           knew,                                  &
     &                           angler, h, pm, pn, on_u, om_v,         &
#ifdef MASKING
     &                           vmask,                                 &
#endif
     &                           zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: knew
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: angle, cff, fac, major, minor, omega, phase, val
      real(r8) :: ramp
#if defined BERING
      real(r8) :: my_area, my_width
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  2D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
      IF (LBC(ieast,isUbar,ng)%acquire.and.                             &
     &    LBC(ieast,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrT,JendT
          BOUNDARY(ng)%ubar_east(j)=0.0_r8
        END DO
        DO j=JstrP,JendT
          BOUNDARY(ng)%vbar_east(j)=0.0_r8
        END DO
      END IF
      IF (LBC(ieast,isUbar,ng)%acquire.and.                             &
     &    LBC(ieast,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Western_Edge(tile)) THEN
        my_area =0.0_r8
        my_width=0.0_r8
        DO j=JstrP,JendP
          my_area=my_area+0.5_r8*(h(Istr-1,j)+h(Istr,j))*               &
     &                           on_u(Istr,j)
          my_width=my_width+on_u(Istr,j)
        END DO
        fac=0.4d6
        DO j=JstrT,JendT
          BOUNDARY(ng)%ubar_west(j)=fac/my_area
        END DO
        DO j=JstrP,JendT
          BOUNDARY(ng)%vbar_west(j)=0.0_r8
        END DO
      END IF
      IF (LBC(ieast,isUbar,ng)%acquire.and.                             &
     &    LBC(ieast,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Southern_Edge(tile)) THEN
        my_area =0.0_r8
        my_width=0.0_r8
        DO i=IstrP,IendP
          my_area=my_area+0.5_r8*(h(i,Jstr-1)+h(i,Jstr))*               &
     &                           om_v(i,Jstr)
          my_width=my_width+om_v(i,Jstr)
        END DO
        fac=0.4d6
        DO i=IstrP,IendT
          BOUNDARY(ng)%ubar_south(i)=0.0_r8
        END DO
        DO i=IstrT,IendT
          BOUNDARY(ng)%vbar_south(i)=fac/my_area
        END DO
      END IF
      IF (LBC(ieast,isUbar,ng)%acquire.and.                             &
     &    LBC(ieast,isVbar,ng)%acquire.and.                             &
     &    DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrP,IendT
          BOUNDARY(ng)%ubar_north(i)=0.0_r8
        END DO
        DO i=IstrT,IendT
          BOUNDARY(ng)%vbar_north(i)=0.0_r8
        END DO
      END IF
      RETURN
      END SUBROUTINE ana_m2obc_tile
