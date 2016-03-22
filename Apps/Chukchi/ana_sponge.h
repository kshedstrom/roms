      SUBROUTINE ana_sponge (ng, tile, model)
!
!! svn $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!                                                                      !
!  This routine rescales horizontal mixing coefficients according      !
!  to the grid size.  Also,  if applicable,  increases horizontal      !
!  in sponge areas.                                                    !
!                                                                      !
!  WARNING:   All biharmonic coefficients are assumed to have the      !
!             square root taken and have  m^2 s^-1/2 units.  This      !
!             will allow multiplying the  biharmonic  coefficient      !
!             to harmonic operator.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"

      CALL ana_sponge_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(8)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_sponge
!
!***********************************************************************
      SUBROUTINE ana_sponge_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_scalars
!
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
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
      integer :: Iwrk, i, ifoo, itwo, j, itrc
      real(r8) :: cff, cff1, cff2, fac

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing in the sponge areas.
!-----------------------------------------------------------------------
!
!! User modifiable section.  Please specify the appropiate sponge area
!! by increasing its horizontal mixing coefficients.
!!            

!
! Chukchi sponge areas
!
      Iwrk = 10
#if defined UV_VIS2
      DO i=IstrT,IendT
        DO j=JstrT,MIN(Iwrk,JendT)
          cff = 250.*0.5_r8*(1.0_r8+COS(pi*REAL(j,r8)/REAL(Iwrk,r8)))
          MIXING(ng) % visc2_r(i,j) = max(cff,                          &
     &                                MIXING(ng) % visc2_r(i,j))
          MIXING(ng) % visc2_p(i,j) = max(cff,                          &
     &                                MIXING(ng) % visc2_p(i,j))
        END DO
      END DO
      DO i=IstrT,MIN(Iwrk,IendT)
        DO j=MAX(JstrT,i),JendT
          cff = 250.*0.5_r8*(1.0_r8+COS(pi*REAL(i,r8)/REAL(Iwrk,r8)))
          MIXING(ng) % visc2_r(i,j) = max(cff,                          &
     &                                MIXING(ng) % visc2_r(i,j))
          MIXING(ng) % visc2_p(i,j) = max(cff,                          &
     &                                MIXING(ng) % visc2_p(i,j))
        END DO
      END DO
!     DO i=MAX(Lm(ng)+1-Iwrk,IstrT),IendT
!       ifoo = Lm(ng)+1-i
!       DO j=MAX(JstrT,ifoo),JendT
!         cff = 250.*0.5_r8*(1.0_r8+COS(pi*REAL(ifoo,r8)/REAL(Iwrk,r8)))
!         MIXING(ng) % visc2_r(i,j) = max(cff,                          &
!    &                                MIXING(ng) % visc2_r(i,j))
!         MIXING(ng) % visc2_p(i+1,j) = max(cff,                        &
!    &                                MIXING(ng) % visc2_p(i+1,j))
!       END DO
!     END DO
      DO j=MAX(Mm(ng)+1-Iwrk,JstrT),JendT
        ifoo = Mm(ng)+1-j
        itwo = Lm(ng)-Mm(ng)+j
        DO i=MAX(IstrT,ifoo),MIN(IendT,itwo)
          cff = 250.*0.5_r8*(1.0_r8+COS(pi*REAL(ifoo,r8)/REAL(Iwrk,r8)))
          MIXING(ng) % visc2_r(i,j) = max(cff,                          &
     &                                MIXING(ng) % visc2_r(i,j))
          MIXING(ng) % visc2_p(i+1,j) = max(cff,                        &
     &                                MIXING(ng) % visc2_p(i+1,j))
        END DO
      END DO
#endif
#ifdef SOLVE3D
# if defined TS_DIF2
      DO itrc=1,NT(ng)
        DO j=JstrT,MIN(Iwrk,JendT)
          cff = 100. * (1.0_r8+COS(pi*REAL(j,r8)/REAL(Iwrk,r8)))
          DO i=IstrT,IendT
            MIXING(ng) % diff2(i,j,itrc)=max(cff,                       &
     &                                MIXING(ng) % diff2(i,j,itrc))
          END DO
        END DO
        DO i=IstrT,MIN(Iwrk,IendT)
          DO j=MAX(JstrT,i),JendT
            cff = 100. * (1.0_r8+COS(pi*REAL(i,r8)/REAL(Iwrk,r8)))
            MIXING(ng) % diff2(i,j,itrc) = max(cff,                     &
     &                                MIXING(ng) % diff2(i,j,itrc))
          END DO
        END DO
!       DO i=MAX(Lm(ng)+1-Iwrk,IstrT),IendT
!         ifoo = Lm(ng)+1-i
!         DO j=MAX(JstrT,ifoo),JendT
!           cff = 100. * (1.0_r8+COS(pi*REAL(ifoo,r8)/REAL(Iwrk,r8)))
!           MIXING(ng) % diff2(i,j,itrc) = max(cff,                     &
!    &                                MIXING(ng) % diff2(i,j,itrc))
!         END DO
!       END DO
        DO j=MAX(Mm(ng)+1-Iwrk,JstrT),JendT
          ifoo = Mm(ng)+1-j
          itwo = Lm(ng)-Mm(ng)+j
          DO i=MAX(IstrT,ifoo),MIN(IendT,itwo)
            cff = 100. * (1.0_r8+COS(pi*REAL(ifoo,r8)/REAL(Iwrk,r8)))
            MIXING(ng) % diff2(i,j,itrc) = max(cff,                     &
     &                                MIXING(ng) % diff2(i,j,itrc))
          END DO
        END DO
      END DO
# endif
#endif
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
!! WARNING:  This section is generic for all applications. Please do not
!!           change the code below.
!!            
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
#ifdef UV_VIS2
        CALL exchange_r2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        MIXING(ng) % visc2_r)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        MIXING(ng) % visc2_p)
#endif
#ifdef UV_VIS4
        CALL exchange_r2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        MIXING(ng) % visc4_r)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        MIXING(ng) % visc4_p)
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
        DO itrc=1,NT(ng)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          MIXING(ng) % diff2(:,:,itrc))
        END DO
# endif
# ifdef TS_DIF4
        DO itrc=1,NT(ng)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          MIXING(ng) % diff4(:,:,itrc))
        END DO
# endif
#endif
      END IF
#ifdef DISTRIBUTE
# ifdef UV_VIS2
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    MIXING(ng) % visc2_r, MIXING(ng) % visc2_p)
# endif
# ifdef UV_VIS4
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    MIXING(ng) % visc4_r, MIXING(ng) % visc4_p)
# endif
# ifdef SOLVE3D
#  ifdef TS_DIF2
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, NT(ng),                &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    MIXING(ng) % diff2)
#  endif
#  ifdef TS_DIF4
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, NT(ng),                &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    MIXING(ng) % diff4)
#  endif
# endif
#endif
      RETURN
      END SUBROUTINE ana_sponge_tile
