      SUBROUTINE ana_snow (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets surface snowfall rate.                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
      USE mod_ice
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_snow_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    ICE(ng) % tis,                                &
     &                    FORCES(ng) % snow)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(17)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_snow
!
!***********************************************************************
      SUBROUTINE ana_snow_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          tis, snow)
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
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: tis(LBi:,LBj:)
      real(r8), intent(out) :: snow(LBi:,LBj:)
#else
      real(r8), intent(in)  :: tis(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: snow(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: iday, month, year
      real(r8) :: hour, yday
      real(r8), parameter :: sn(12) =                                   &
     &        (/ 1.1e-6, 1.1e-6, 1.1e-6, 1.1e-6, 6.4e-6, 0.,            &
     &           0., 0., 1.64e-5, 1.64e-5, 1.1e-6, 1.1e-6    /)
!
! Snow from Maykut and Untersteiner:
!   30 cm between Aug 20 and Oct 30, 5 cm to April 30  + 5 cm in May.
! Assuming dry snow density of 330 kg/m^3
!

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set snow precipitation rate (kg/m2/s).
!-----------------------------------------------------------------------
!
#ifdef NO_SNOW
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          snow(i,j)=0.0_r8
        END DO
      END DO
#else
      CALL caldate(r_date, tdays(ng), year, yday, month, iday, hour)
      IF (month == 8 .and. iday >= 20) month = 9
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          snow(i,j)=sn(month)
!          snow(i,j)=sn(month) * 2.0_r8
          IF (tis(i,j) > 0.0_r8) snow(i,j) = 0.0_r8
        END DO
      END DO
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          snow)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    snow)
#endif

      RETURN
      END SUBROUTINE ana_snow_tile
