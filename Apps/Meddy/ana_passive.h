      SUBROUTINE ana_passive (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets initial conditions for passive inert tracers      !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_ocean
      USE mod_grid
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_passive_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
#ifdef SPHERICAL
     &                       GRID(ng) % lonr,                           &
     &                       GRID(ng) % latr,                           &
#else
     &                       GRID(ng) % xr,                             &
     &                       GRID(ng) % yr,                             &
#endif
     &                       GRID(ng) % z_r,                            &
     &                       OCEAN(ng) % t)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(18)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_passive
!
!***********************************************************************
      SUBROUTINE ana_passive_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
#ifdef SPHERICAL
     &                             lonr, latr,                          &
#else
     &                             xr, yr,                              &
#endif
     &                             z_r,                                 &
     &                             t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
# else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
# endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, ip, itrc, j, k
      real(r8) :: val1

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical initial conditions for passive inert tracers.
!-----------------------------------------------------------------------
!
# if defined MEDDY
      DO ip=1,NPT
        itrc=inert(ip)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              t(i,j,k,1,itrc) = 0.
            END DO
          END DO
        END DO
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              val1 = sqrt( ((xr(i,j)-40.e+3_r8)/10000.)**2 +              &
     &                     ((yr(i,j)-40.e+3_r8)/10000.)**2 +              &
     &                     ((z_r(i,j,k) + 500.)/200.)**2 )
              IF (val1 .le. 1) THEN
                t(i,j,k,1,itrc)=1.0
              END IF
              t(i,j,k,2,itrc)=t(i,j,k,1,itrc)
            END DO
          END DO
        END DO
      END DO
#else
      ana_passive_user.h: No values provided for passive tracers.
#endif

      RETURN
      END SUBROUTINE ana_passive_tile
