      SUBROUTINE ana_tclima (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets analytical tracer climatology fields.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_tclima_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      CLIMA(ng) % tclm)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(33)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_tclima
!
!***********************************************************************
      SUBROUTINE ana_tclima_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            tclm)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
      USE mod_boundary
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange4d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: tclm(LBi:,LBj:,:,:)
#else
      real(r8), intent(out) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k, ifoo

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set tracer climatology.
!-----------------------------------------------------------------------
!
#if defined CCS1
      ifoo = 10
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            tclm(i,j,k,itemp)=0.0
            tclm(i,j,k,isalt)=0.0
          END DO
        END DO
        DO j=JstrT,MIN(ifoo,JendT)               ! SOUTH boundary
          DO i=IstrT,IendT
            tclm(i,j,k,itemp)=BOUNDARY(ng)%t_south(i,k,itemp)
            tclm(i,j,k,isalt)=BOUNDARY(ng)%t_south(i,k,isalt)
          END DO 
        END DO
        DO j=MAX(JstrT,Mm(ng)+1-ifoo),JendT       ! NORTH boundary
          DO i=IstrT,IendT
            tclm(i,j,k,itemp)=BOUNDARY(ng)%t_north(i,k,itemp)
            tclm(i,j,k,isalt)=BOUNDARY(ng)%t_north(i,k,isalt)
          END DO
        END DO
        DO i=IstrT,MIN(ifoo,IendT)                ! WEST boundary
          DO j=JstrT,JendT
            tclm(i,j,k,itemp)=BOUNDARY(ng)%t_west(j,k,itemp)
            tclm(i,j,k,isalt)=BOUNDARY(ng)%t_west(j,k,isalt)
          END DO
        END DO
      END DO
#else
      ana_tclima.h: No values provided for tclm.
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tclm(:,:,:,itrc))
        END DO
      END IF
#ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tclm)
#endif

      RETURN
      END SUBROUTINE ana_tclima_tile
