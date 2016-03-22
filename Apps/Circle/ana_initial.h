      SUBROUTINE ana_initial (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This subroutine sets initial conditions for momentum and tracer     !
!  type variables using analytical expressions.                        !
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
      CALL ana_initial_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % xr,                             &
     &                       GRID(ng) % yr,                             &
     &                       GRID(ng) % xu,                             &
     &                       GRID(ng) % yu,                             &
     &                       GRID(ng) % xv,                             &
     &                       GRID(ng) % yv,                             &
#ifdef MASKING
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
     &                       GRID(ng) % rmask,                          &
#endif
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(10)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_initial
!
!***********************************************************************
      SUBROUTINE ana_initial_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             h,                                   &
     &                             xr, yr,                              &
     &                             xu, yu, xv, yv,                      &
#ifdef MASKING
     &                             umask, vmask, rmask,                 &
#endif
     &                             ubar, vbar, zeta)
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
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
      real(r8), intent(in) :: xu(LBi:,LBj:)
      real(r8), intent(in) :: yu(LBi:,LBj:)
      real(r8), intent(in) :: xv(LBi:,LBj:)
      real(r8), intent(in) :: yv(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(out) :: ubar(LBi:,LBj:,:)
      real(r8), intent(out) :: vbar(LBi:,LBj:,:)
      real(r8), intent(out) :: zeta(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: xu(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yu(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: xv(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yv(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: Iless, Iplus, i, itrc, j, k, nu, NM
      real(r8) :: rad, f0, sigma, omega, grav, cu, bessi_c
      real(r8) :: theta, rad_tilde, Fx, Fy, kappa, Fxc, Fxs, Fyd
      real(r8) :: FF, vt, vr, depth, c0

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initial conditions for 2D momentum (m/s) components.
!-----------------------------------------------------------------------
!
! Some constants
      grav = 3.92e-2
      f0 = 8.34e-5
      depth = 1500.0
      c0 = 1.                      ! wave amplitude
      nu = 3
      omega = 0.4986
      sigma = omega*f0
      kappa = sqrt( (f0**2 - sigma**2)/(grav*depth) )

#if defined CIRCLE
      DO j=JstrT,JendT
        DO i=IstrP,IendT
! u-velocity points
	  IF (umask(i,j) .gt. 0.5) THEN
            theta = atan2(yu(i,j),xu(i,j))
            rad = sqrt(xu(i,j)**2 + yu(i,j)**2)
            rad_tilde = rad*kappa
            Fxc = cos(theta*nu)
            Fxs = sin(theta*nu)
            Fy = bessi_c(nu,rad_tilde)
            Fyd = bessi_c(nu-1,rad_tilde)

            FF = kappa*Fyd - nu/rad * Fy
            vt = c0*Fxc/(kappa**2*depth) * (f0*FF - sigma*nu*Fy/rad)
            vr = -c0*Fxs/(kappa**2*depth) * (sigma*FF - nu*f0*Fy/rad)
   
            ubar(i,j,1) = vr*cos(theta) - vt*sin(theta)
          ELSE
            ubar(i,j,1) = 0.0
          END IF
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
! v-velocity points
	  IF (vmask(i,j) .gt. 0.5) THEN
            theta = atan2(yv(i,j),xv(i,j))
            rad = sqrt(xv(i,j)**2 + yv(i,j)**2)
            rad_tilde = rad*kappa
            Fxc = cos(theta*nu)
            Fxs = sin(theta*nu)
            Fy = bessi_c(nu,rad_tilde)
            Fyd = bessi_c(nu-1,rad_tilde)

            FF = kappa*Fyd - nu/rad * Fy
            vt = c0*Fxc/(kappa**2*depth) * (f0*FF - sigma*nu*Fy/rad)
            vr = -c0*Fxs/(kappa**2*depth) * (sigma*FF - nu*f0*Fy/rad)

            vbar(i,j,1) = vr*sin(theta) + vt*cos(theta)
          ELSE
            vbar(i,j,1) = 0.0
          END IF
        END DO
      END DO
#else
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ubar(i,j,1)=0.0_r8
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          vbar(i,j,1)=0.0_r8
        END DO
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Initial conditions for free-surface (m).
!-----------------------------------------------------------------------
!
#if defined CIRCLE
      DO j=JstrT,JendT
        DO i=IstrT,IendT
	  IF (rmask(i,j) .gt. 0.5) THEN
            theta = atan2(yr(i,j),xr(i,j))
            rad = sqrt(xr(i,j)**2 + yr(i,j)**2)
            rad_tilde = rad*kappa
            Fx = cos(theta*nu)
            Fy =  bessi_c(nu,rad_tilde)
            zeta(i,j,1) = c0*Fx*Fy
          ELSE
            zeta(i,j,1) = 0.0
          END IF
        END DO
      END DO
#else
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          zeta(i,j,1)=0.0_r8
        END DO
      END DO
#endif

      RETURN
      END SUBROUTINE ana_initial_tile
