      SUBROUTINE ana_stflux (ng, tile, model, itrc)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets kinematic surface flux of tracer type variables   !
!  "stflx" (tracer units m/s) using analytical expressions.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc

#include "tile.h"
!
      CALL ana_stflux_tile (ng, tile, model, itrc,                      &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
#ifdef SHORTWAVE
     &                      FORCES(ng) % srflx,                         &
#endif
#ifdef TL_IOMS
     &                      FORCES(ng) % tl_stflx,                      &
#endif
     &                      FORCES(ng) % stflx)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(31)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_stflux
!
!***********************************************************************
      SUBROUTINE ana_stflux_tile (ng, tile, model, itrc,                &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
#ifdef SHORTWAVE
     &                            srflx,                                &
#endif
#ifdef TL_IOMS
     &                            tl_stflx,                             &
#endif
     &                            stflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
# ifdef SHORTWAVE
      real(r8), intent(in) :: srflx(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
# ifdef TL_IOMS
      real(r8), intent(inout) :: tl_stflx(LBi:,LBj:,:)
# endif
#else
# ifdef SHORTWAVE
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
# ifdef TL_IOMS
      real(r8), intent(inout) :: tl_stflx(LBi:UBi,LBj:UBj,NT(ng))
# endif
#endif
!
!  Local variable declarations.
!
      integer :: i, j
#if defined GAK1D
      integer :: iday, month, year
      real(r8) :: cff, hour, yday
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic surface heat flux (degC m/s) at horizontal
!  RHO-points.
!-----------------------------------------------------------------------
!
      IF (itrc.eq.itemp) THEN
#if defined GAK1D
!       Eyeball fit to COADS Climatological net heat flux near GAK1
        CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
        cff = -125.0_r8 * COS( (yday-24.5_r8) * 2.4_r8  * pi / 360._r8) &
     &         / (rho0*Cp)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            stflx(i,j,itrc)=cff
          END DO
        END DO
#else
        DO j=JstrT,JendT
          DO i=IstrT,IendT
# ifdef BL_TEST
            stflx(i,j,itrc)=srflx(i,j)
#  ifdef TL_IOMS
            tl_stflx(i,j,itrc)=srflx(i,j)
#  endif
# else
            stflx(i,j,itrc)=0.0_r8
#  ifdef TL_IOMS
            tl_stflx(i,j,itrc)=0.0_r8
#  endif
# endif
          END DO
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Set kinematic surface freshwater flux (m/s) at horizontal
!  RHO-points, scaling by surface salinity is done in STEP3D.
!-----------------------------------------------------------------------
!
      ELSE IF (itrc.eq.isalt) THEN
#ifdef GAK1D
!       Tuned to generate S profile at GAK1 - includes effect of runoff
        CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
        cff = ( -0.025_r8 + 0.025_r8 * &
     &    COS( (yday-61._r8) * 2.0_r8 * pi / 360._r8 ) ) / 86400.0_r8
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            stflx(i,j,itrc) = cff
          END DO
        END DO
#else
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            stflx(i,j,itrc)=0.0_r8
# ifdef TL_IOMS
            tl_stflx(i,j,itrc)=0.0_r8
# endif
          END DO
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Set kinematic surface flux (T m/s) of passive tracers, if any.
!-----------------------------------------------------------------------
!
      ELSE
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            stflx(i,j,itrc)=0.0_r8
#ifdef TL_IOMS
            tl_stflx(i,j,itrc)=0.0_r8
#endif
          END DO
        END DO
      END IF
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        stflx(:,:,itrc))
#ifdef TL_IOMS
        CALL exchange_r2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        tl_stflx(:,:,itrc))
#endif
      END IF
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    stflx(:,:,itrc))
# ifdef TL_IOMS
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_stflx(:,:,itrc))
# endif
#endif
      RETURN
      END SUBROUTINE ana_stflux_tile
