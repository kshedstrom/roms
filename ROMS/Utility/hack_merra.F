#include "cppdefs.h"
      Module hack_mod

      Contains
#ifdef ARCTIC_MERRA_HACK
!***********************************************************************
      SUBROUTINE hack_flux_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           srflx, lrflx, ai)
!***********************************************************************
      USE mod_param
      USE mod_scalars
      USE mod_stepping
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj 
# ifdef ASSUMED_SHAPE
#  ifdef CICE_MODEL
      real(r8), intent(in) :: ai(LBi:,LBj:)
#  else
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
# else
#  ifdef CICE_MODEL
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
#  endif
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lrflx(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: iday, month, year
      real(r8) :: hour, yday

!----------------------------------------------------------------------
! Values for Barrow from Walsh et al., 2009, J. Climate.
!
! Cullather and Bosilovich, 2012, J. Climate would say to add 10 more
! each month to longwave and to cool spring Tair to match SHEBA.
!----------------------------------------------------------------------
!      real(r8), parameter :: lw_off(12) =                               &
!     &        (/ -5, -10, -22, -15, 14, 10, 10, 12, 20, 20, -2, -5 /)
! SHEBA correction to the Walsh version
      real(r8), parameter :: lw_off(12) =                               &
     &        (/ 5, 0, -12, -5, 24, 20, 20, 22, 30, 30, 8, 5 /)
      real(r8), parameter :: sw_off(12) =                               &
     &        (/ 0, 10, 25, 50, 25, 0, -5, -10, -5, 0, 0, 0 /)
      real(r8) :: fac1, fac2
!
!----------------------------------------------------------------------
!
# include "set_bounds.h"
!
!----------------------------------------------------------------------
!
      CALL caldate(r_date, tdays(ng), year, yday, month, iday, hour)
      fac1 = lw_off(month)/(rho0*Cp)
      fac2 = sw_off(month)/(rho0*Cp)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
#  ifdef CICE_MODEL
          lrflx(i,j)=lrflx(i,j) + ai(i,j)*fac1
          srflx(i,j)=MAX(0.0_r8, srflx(i,j) + ai(i,j)*fac2)
#  else
          lrflx(i,j)=lrflx(i,j) + ai(i,j,linew(ng))*fac1
          srflx(i,j)=MAX(0.0_r8, srflx(i,j) + ai(i,j,linew(ng))*fac2)
#  endif
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          srflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          lrflx)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    srflx, lrflx)
# endif
      RETURN
      END SUBROUTINE hack_flux_tile
#endif
      end module
