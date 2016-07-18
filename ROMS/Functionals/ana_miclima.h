      SUBROUTINE ana_miclima (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets analytical ice momentum climatology fields.       !
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
      CALL ana_miclima_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(11)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_miclima
!
!***********************************************************************
      SUBROUTINE ana_miclima_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
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
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set 2D momentum climatology.
!-----------------------------------------------------------------------
!
      IF (LmiCLM(ng)) THEN
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            CLIMA(ng)%uiclm(i,j)=???
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            CLIMA(ng)%viclm(i,j)=???
          END DO
        END DO
!
!  Exchange boundary data.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            CLIMA(ng) % uiclm)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            CLIMA(ng) % viclm)
        END IF

#ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      CLIMA(ng) % uiclm,                          &
     %                      CLIMA(ng) % viclm)
#endif
      END IF

      RETURN
      END SUBROUTINE ana_miclima_tile
