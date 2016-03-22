      SUBROUTINE ana_m3obc (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets 3D momentum open boundary conditions using        !
!  analytical expressions.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_m3obc_tile (ng, tile, model,                             &
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
        ANANAME(14)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_m3obc
!
!***********************************************************************
      SUBROUTINE ana_m3obc_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: fac, val

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  3D momentum open boundary conditions.
!-----------------------------------------------------------------------
!
#if defined SED_TEST1
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        fac=5.0E-06_r8
        DO k=1,N(ng)
          DO j=JstrT,JendT
            val=0.5_r8*(zeta(0 ,j,knew)+h(0 ,j)+                        &
     &                  zeta(1 ,j,knew)+h(1 ,j))
            BOUNDARY(ng)%u_west(j,k)=-LOG((val+0.5*(z_r(Istr-1,j,k)+    &
     &                                              z_r(Istr  ,j,k)))/  &
     &                                    fac)/                         &
     &                               (LOG(val/fac)-1.0_r8+fac/val)
          END DO
          DO j=Jstr,JendT
            BOUNDARY(ng)%v_west(j,k)=0.0_r8
          END DO
        END DO
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        fac=5.0E-06_r8
        DO k=1,N(ng)
          DO j=JstrT,JendT
            val=0.5_r8*(zeta(Iend  ,j,knew)+h(Iend  ,j)+                &
     &                  zeta(Iend+1,j,knew)+h(Iend+1,j))
            BOUNDARY(ng)%u_east(j,k)=-LOG((val+0.5*(z_r(Iend  ,j,k)+    &
     &                                              z_r(Iend+1,j,k)))/  &
     &                                    fac)/                         &
     &                               (LOG(val/fac)-1.0_r8+fac/val)
          END DO
          DO j=Jstr,JendT
            BOUNDARY(ng)%v_east(j,k)=0.0_r8
          END DO
        END DO
      END IF
#elif defined LA_04
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
       DO k=1,N(ng)
        DO j=JstrT,JendT
	  BOUNDARY(ng)%u_west(j,k)=0.0179588894377609_r8
        END DO
	DO j=Jstr,JendT
	  BOUNDARY(ng)%v_west(j,k)=0.0_r8
        END DO
       END DO
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO k=1,N(ng)
         DO j=JstrT,JendT
	   BOUNDARY(ng)%u_east(j,k)=3.72297256457142E-5_r8
         END DO
	 DO j=JstrP,JendT
	   BOUNDARY(ng)%v_east(j,k)=0.0_r8
         END DO
        END DO
      END IF
#else
# ifdef EAST_M3OBS
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            BOUNDARY(ng)%u_east(j,k)=0.0_r8
          END DO
          DO j=JstrP,JendT
            BOUNDARY(ng)%v_east(j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
# ifdef WEST_M3OBS
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            BOUNDARY(ng)%u_west(j,k)=0.0_r8
          END DO
          DO j=Jstr,JendT
            BOUNDARY(ng)%v_west(j,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
# ifdef SOUTH_M3OBS
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=IstrP,IendT
            BOUNDARY(ng)%u_south(i,k)=0.0_r8
          END DO
          DO i=IstrT,IendT
            BOUNDARY(ng)%v_south(i,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
# ifdef NORTH_M3OBS
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO k=1,N(ng)
          DO i=IstrP,IendT
            BOUNDARY(ng)%u_north(i,k)=0.0_r8
          END DO
          DO i=IstrT,IendT
            BOUNDARY(ng)%v_north(i,k)=0.0_r8
          END DO
        END DO
      END IF
# endif
#endif
      RETURN
      END SUBROUTINE ana_m3obc_tile
