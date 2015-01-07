      SUBROUTINE ana_nudgcoef (ng, tile, model)
!
!! svn $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets spatially varying nudging coefficients time-      !
!  scales (1/s). They are used for nuding to climatology in the        !
!  governing equations.                                                !
!                                                                      !
!  It is HIGHLY recommended to write these nudging coefficients into   !
!  input NetCDF NUDNAME instead of using analytical expressions        !
!  below.  It is very easy to introduce parallel bugs.  Also, Users    !
!  can plot their spatial distribution and fine tune their values      !
!  during the pre-proccessing stage for a particular application.      !
!                                                                      !
!  REMARK:  Nudging of free-surface in the vertically integrated       !
!  ======   continuity equation is NOT allowed because it VIOLATES     !
!  mass/volume conservation. If such nudging effects are required,     !
!  it needs to be specified on the momentum equations for (u,v)        !
!  and/or (ubar,vbar). If done on (u,v) only, its effects enter        !
!  the 2D momentum equations via the residual vertically integrated    !
!  forcing term.                                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL ana_nudgcoef_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(16)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_nudgcoef
!
!***********************************************************************
      SUBROUTINE ana_nudgcoef_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_clima
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_collect
      USE mp_exchange_mod, ONLY : mp_exchange2d
# ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE mp_exchange_mod, ONLY : mp_exchange4d
# endif
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Iwrk, i, itrc, j, k

      real(r8) :: cff1, cff2, cff3

      real(r8), parameter :: IniVal = 0.0_r8

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set up nudging towards data time-scale coefficients (1/s).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          wrk(i,j)=0.0_r8
        END DO
      END DO

!
!  Set nudging boundaries coefficients zone for ARCTIC2 
!  nudging coefficients vary from a thirty
!  days time scale at the boundary point to decrease linearly to 0 days
!  (i.e no nudging) 15 grids points away from the boundary.
!
      cff1=1.0_r8/(30_r8*86400.0_r8)         ! 30-day outer limit
      cff2=0.0_r8                            ! Inf days (no nudge) inner limit
      cff3=20.0_r8                           ! width of layer in grid points

! cff3-point wide linearly tapered nudging zone
      DO j=JstrT,MIN(INT(cff3),JendT)               ! SOUTH boundary
        DO i=IstrT,IendT
          wrk(i,j)=cff2+(cff3-REAL(j,r8))*(cff1-cff2)/cff3
        END DO
      END DO
! cff3-point wide linearly tapered nudging zone
      DO j=MAX(JstrT,Mm(ng)+1-INT(cff3)),JendT       ! NORTH boundary
        DO i=IstrT,IendT
          wrk(i,j)=MAX(wrk(i,j),                                        &
     &             cff1+REAL(Mm(ng)+1-j,r8)*(cff2-cff1)/cff3)
        END DO
      END DO
! cff3-point wide linearly tapered nudging zone
      DO i=IstrT,MIN(INT(cff3),IendT)                ! WEST boundary
        DO j=JstrT,MIN(JendT,120)
          wrk(i,j)=MAX(wrk(i,j),                                        &
     &             cff2+(cff3-REAL(i,r8))*(cff1-cff2)/cff3)
        END DO
      END DO
! cff3-point wide linearly tapered nudging zone
      DO i=MAX(IstrT,Lm(ng)+1-INT(cff3)),IendT       ! EAST boundary
        DO j=JstrT,MIN(JendT,410)
          wrk(i,j)=MAX(wrk(i,j),                                        &
     &             cff1+REAL(Lm(ng)+1-i,r8)*(cff2-cff1)/cff3)
        END DO
      END DO
!
! Set the relevant nudging coefficients using the entries in wrk
!
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        DO itrc=1,NTCLM(ng)
          DO k=1,N(ng)
            DO j=JstrT,JendT
              DO i=IstrT,IendT
                CLIMA(ng)%Tnudgcof(i,j,k,itrc)=wrk(i,j)
              END DO
            END DO
          END DO
        END DO
      END IF
      IF (LnudgeM2CLM(ng)) THEN
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            CLIMA(ng)%M2nudgcof(i,j)=wrk(i,j)
          END DO
        END DO
      END IF
      IF (LnudgeM3CLM(ng)) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              CLIMA(ng)%M3nudgcof(i,j,k)=wrk(i,j)
            END DO
          END DO
        END DO
      END IF

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchage nudging coefficients information.
!-----------------------------------------------------------------------
!
      IF (LnudgeM2CLM(ng)) THEN
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      CLIMA(ng)%M2nudgcof)
      END IF

# ifdef SOLVE3D
!
      IF (LnudgeM3CLM(ng)) THEN
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      CLIMA(ng)%M3nudgcof)
      END IF
!
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        CALL mp_exchange4d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng), 1, NTCLM(ng), &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      CLIMA(ng)%Tnudgcof)
      END IF
# endif
#endif

      RETURN
      END SUBROUTINE ana_nudgcoef_tile
