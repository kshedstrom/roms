     SUBROUTINE ice_flux_rst (ng, tile)
!
!*************************************************** W. Paul Budgell ***
!  Copyright (c) 2002-2010 ROMS/TOMS Group                             !
!************************************************** Hernan G. Arango ***
!                                                                      !
!  The purpose of this routine is to provide the ice-water fluxes      !
!  during the first step of a perfect restart.                         !
!                                                                      !
!***********************************************************************
!

      USE mod_param
      USE mod_ice
      USE mod_forces

      implicit none

      integer, intent(in) :: ng, tile

#include "tile.h"

#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 51)
#endif
      CALL ice_flux_rst_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        FORCES(ng) % sustr_save,                  &
     &                        FORCES(ng) % svstr_save,                  &
     &                        FORCES(ng) % stflx_save,                  &
     &                        FORCES(ng) % sustr,                       &
     &                        FORCES(ng) % svstr,                       &
     &                        FORCES(ng) % stflx)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 51)
#endif
      RETURN
      END SUBROUTINE ice_flux_rst
!
!***********************************************************************
      SUBROUTINE ice_flux_rst_tile (ng, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        sustr_save, svstr_save,                   &
     &                        stflx_save, sustr, svstr,                 &
     &                        stflx)
!***********************************************************************
!
!***********************************************************************

      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: sustr_save(LBi:,LBj:,:)
      real(r8), intent(in) :: svstr_save(LBi:,LBj:,:)
      real(r8), intent(in) :: stflx_save(LBi:,LBj:,:)
      real(r8), intent(out) :: sustr(LBi:,LBj:,:)
      real(r8), intent(out) :: svstr(LBi:,LBj:,:)
      real(r8), intent(out) :: stflx(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: sustr_save(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: svstr_save(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: stflx_save(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(out) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
#endif

! Local variable definitions
!
      integer :: i, j, itrc

#include "set_bounds.h"

      DO itrc=1,NAT
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            sustr(i,j,itrc) = sustr_save(i,j,itrc)
            svstr(i,j,itrc) = svstr_save(i,j,itrc)
            stflx(i,j,itrc) = stflx_save(i,j,itrc)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE ice_flux_rst_tile
