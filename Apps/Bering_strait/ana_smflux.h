      SUBROUTINE ana_smflux (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets kinematic surface momentum flux (wind stress)     !
!  "sustr" and "svstr" (m2/s2) using an analytical expression.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_smflux_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng) % angler,                          &
     &                      GRID(ng) % mask2,                           &
#ifdef SPHERICAL
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
#else
     &                      GRID(ng) % xr,                              &
     &                      GRID(ng) % yr,                              &
#endif
#ifdef TL_IOMS
     &                      FORCES(ng) % tl_sustr,                      &
     &                      FORCES(ng) % tl_svstr,                      &
#endif
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(24)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_smflux
!
!***********************************************************************
      SUBROUTINE ana_smflux_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            angler,                               &
     &                            mask2,                                &
#ifdef SPHERICAL
     &                            lonr, latr,                           &
#else
     &                            xr, yr,                               &
#endif
#ifdef TL_IOMS
     &                            tl_sustr, tl_svstr,                   &
#endif
     &                            sustr, svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: mask2(LBi:,LBj:)
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
# else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
# endif
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
# ifdef TL_IOMS
      real(r8), intent(out) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(out) :: tl_svstr(LBi:,LBj:)
# endif
#else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: mask2(LBi:UBi,LBj:UBj)
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
# ifdef TL_IOMS
      real(r8), intent(out) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_svstr(LBi:UBi,LBj:UBj)
# endif
#endif
!
!  Local variable declarations.
!
      integer :: i, j, i0, j0
      real(r8) :: Ewind, Nwind, windamp, winddir, radius, rad0
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: speed

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (wind stress) component in the
!  XI-direction (m2/s2) at horizontal U-points.
!-----------------------------------------------------------------------
!
#define GYRE
#ifdef GYRE
      i0 = 293
      j0 = 97
      rad0 = 100._r8
      windamp = 0.0001_r8
      DO j=JstrT,JendT
        DO i=IstrT,IendT
! Weighting j variations more
          radius = sqrt((i-i0)**2 + 1._r8*(j-j0)**2)
          speed(i,j) = windamp*0.5*radius / rad0 *                      &
     &                  (tanh((radius-rad0)/rad0) - 1.0)
        END DO
      END DO
      DO j=JstrT,JendT
        DO i=IstrT,IendT
	  winddir = atan2(1._r8*(j-j0),1._r8*(i-i0))
          sustr(i,j)= speed(i,j)*sin(winddir)
        END DO
      END DO
      DO j=JstrT,JendT
        DO i=IstrT,IendT
	  winddir = atan2(1._r8*(j-j0),1._r8*(i-i0))
          svstr(i,j)= -speed(i,j)*cos(winddir)
        END DO
      END DO
#elif defined UP_DOWN
! Something that ramps up, then down.
      windamp = 1.0e-4_r8*0.5*(tanh((time(ng) - 86400._r8)/43200._r8)   &
     &                      - tanh((time(ng) - 4*86400._r8)/43200._r8) )
!! This is in degrees clockwise from north
      winddir = 120*pi/180._r8
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          i0 = max(i-1,IstrT)
          sustr(i,j)= windamp*sin(winddir + angler(i,j)) *              &
     &                max(mask2(i,j),mask2(i0,j))
        END DO
      END DO
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          j0 = max(j-1,JstrT)
          svstr(i,j)= windamp*cos(winddir + angler(i,j)) *              &
     &                max(mask2(i,j),mask2(i,j0))
        END DO
      END DO
#else
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          sustr(i,j)=0.0
        END DO
      END DO
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          svstr(i,j)=0.0
!          j0 = max(j-1,JstrT)
!          svstr(i,j)= windamp*cos(winddir + angler(i,j)) *              &
!     &                max(mask2(i,j),mask2(i,j0))
#ifdef TL_IOMS
          tl_svstr(i,j)=Nwind
#endif
        END DO
      END DO
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
#ifdef TL_IOMS
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_sustr)
#endif
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
#ifdef TL_IOMS
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_svstr)
#endif
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
# ifdef TL_IOMS
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_sustr, tl_svstr)
# endif
#endif

      RETURN
      END SUBROUTINE ana_smflux_tile
