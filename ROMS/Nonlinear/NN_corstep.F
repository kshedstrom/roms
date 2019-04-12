#include "cppdefs.h"
      MODULE NN_corstep_mod
#if defined NONLINEAR && defined NN_MIXING && defined SOLVE3D
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine perfoms the corrector step for turbulent kinetic       !
!  energy and length scale prognostic variables, tke and gls.          !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Mellor, G.L. and T. Yamada, 1982:  Development of turbulence        !
!    closure model for geophysical fluid problems, Rev. Geophys.       !
!    Space Phys., 20, 851-875.                                         !
!                                                                      !
!  Galperin, B., L.H. Kantha, S. Hassid, and A.Rosati, 1988:  A        !
!    quasi-equilibrium  turbulent  energy model for geophysical        !
!    flows, J. Atmos. Sci., 45, 55-62.                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: NN_corstep

      CONTAINS
!
!***********************************************************************
      SUBROUTINE NN_corstep (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 20, __LINE__, __FILE__)
# endif
      CALL NN_corstep_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng), nnew(ng),                         &
# ifdef MASKING
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
     &                      GRID(ng) % Huon,                            &
     &                      GRID(ng) % Hvom,                            &
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
     &                      OCEAN(ng) % W,                              &
# ifdef JK_MIXING
     &                      OCEAN(ng) % rho,                            &
     &                      OCEAN(ng) % pden,                           &
# endif
     &                      FORCES(ng) % bustr,                         &
     &                      FORCES(ng) % bvstr,                         &
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr,                         &
     &                      MIXING(ng) % bvf,                           &
     &                      MIXING(ng) % Akt,                           &
     &                      MIXING(ng) % Akv,                           &
     &                      MIXING(ng) % Akk,                           &
     &                      MIXING(ng) % Lscale,                        &
     &                      MIXING(ng) % gls,                           &
     &                      MIXING(ng) % tke)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 20, __LINE__, __FILE__)
# endif
      RETURN
      END SUBROUTINE NN_corstep
!
!***********************************************************************
      SUBROUTINE NN_corstep_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp, nnew,                           &
# ifdef MASKING
     &                            umask, vmask,                         &
# endif
     &                            Huon, Hvom, Hz, pm, pn, z_r, z_w,     &
     &                            u, v, W,                              &
# ifdef JK_MIXING
     &                            rho, pden,                            &
# endif
     &                            bustr, bvstr, sustr, svstr,           &
     &                            bvf, Akt, Akv,                        &
     &                            Akk, Lscale, gls, tke)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
# ifdef OVETTURNS
      USE rho_eos_mod, ONLY : rho_eos_tile
# endif
      USE tkebc_mod, ONLY : tkebc_tile
      USE akbc_mod, ONLY : akbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
# ifdef JK_MIXING
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: pden(LBi:,LBj:,:)
# endif
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)

      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akk(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Lscale(LBi:,LBj:,0:)
      real(r8), intent(inout) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tke(LBi:,LBj:,0:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
# ifdef JK_MIXING
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pden(LBi:UBi,LBj:UBj,N(ng))
# endif
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))

      real(r8), intent(inout) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(inout) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: Akk(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: Lscale(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: gls(LBi:UBi,LBj:UBj,0:N(ng),3)
      real(r8), intent(inout) :: tke(LBi:UBi,LBj:UBj,0:N(ng),3)
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8), parameter :: Gadv = 1.0_r8/3.0_r8
      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8) :: Gh, Ls_unlmt, Ls_lmt, Qprod, Qdiss, Sh, Sm, Wscale
      real(r8) :: cff, cff1, cff2, cff3, ql, strat2
      real(r8) :: cffm1, cffm2, cffm3, cffm4, cffm5, cffm6
      real(r8) :: alpha_c, Gm, phi1, phi2, phi3, phi4, phi5, Sq
      real(r8) :: q, qint, qzint, lbinv, ls, dum

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BCK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BCP
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCP
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dU
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dV

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: shear2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: buoy2

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FEK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FEP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FXK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FXP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curvK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curvP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradP
      real(r8), dimension(IminS:ImaxS) :: lt
      real(r8), dimension(0:N(ng)) :: AKdum
# ifdef RR_MIXING
      real(r8) :: Riinv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: SmSh
# endif
# ifdef JK_MIXING
      real(r8) :: diffdum, pdenmax
      integer :: kdum, kk
      logical :: found
      real(r8), dimension(0:N(ng)) :: L_t
      real(r8), dimension(0:N(ng)) :: rho_sort
      real(r8), dimension(0:N(ng)) :: bvf_sort
      real(r8), dimension(0:N(ng)) :: rho_dum
      real(r8), dimension(0:N(ng)) :: rho_dums
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: SmSh
# endif

# include "set_bounds.h"
!
!   Compute some constants for NN
       cffm1=-9.0_r8*my_A1*my_A2*(1.0_r8-my_C2)
       cffm2=9.0_r8*my_A2*my_A2*(1.0_r8-my_C2)*(1.0_r8-my_C5)
       cffm3=-12.0_r8*my_A1*my_A2*(1.0_r8-my_C2)
       cffm4=6.0_r8*my_A1*my_A1
       cffm5=-3.0_r8*my_A2*my_B2*(1.0_r8-my_C3)
       cffm6=1.0_r8/sqrt(2.0_r8)
!
!-----------------------------------------------------------------------
!  Compute vertical velocity shear at W-points.
!-----------------------------------------------------------------------
!
# ifdef RI_SPLINES
      DO j=Jstrm1,Jendp1
        DO i=Istrm1,Iendp1
          CF(i,0)=0.0_r8
          dU(i,0)=0.0_r8
          dV(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=Istrm1,Iendp1
            cff=1.0_r8/(2.0_r8*Hz(i,j,k+1)+                             &
     &                  Hz(i,j,k)*(2.0_r8-CF(i,k-1)))
            CF(i,k)=cff*Hz(i,j,k+1)
            dU(i,k)=cff*(3.0_r8*(u(i  ,j,k+1,nstp)-u(i,  j,k,nstp)+     &
     &                           u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))-    &
     &                   Hz(i,j,k)*dU(i,k-1))
            dV(i,k)=cff*(3.0_r8*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+     &
     &                           v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp))-    &
     &                   Hz(i,j,k)*dV(i,k-1))
          END DO
        END DO
        DO i=Istrm1,Iendp1
          dU(i,N(ng))=0.0_r8
          dV(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=Istrm1,Iendp1
            dU(i,k)=dU(i,k)-CF(i,k)*dU(i,k+1)
            dV(i,k)=dV(i,k)-CF(i,k)*dV(i,k+1)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istrm1,Iendp1
            shear2(i,j,k)=dU(i,k)*dU(i,k)+dV(i,k)*dV(i,k)
          END DO
        END DO
      END DO
# else
      DO k=1,N(ng)-1
        DO j=Jstrm1,Jendp1
          DO i=Istrm1,Iendp1
            cff=0.5_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            shear2(i,j,k)=(cff*(u(i  ,j,k+1,nstp)-u(i  ,j,k,nstp)+      &
     &                          u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp)))**2+ &
     &                    (cff*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+      &
     &                          v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)))**2
          END DO
        END DO
      END DO
# endif
!
! Load Brunt-Vaisala frequency.
!
      DO k=1,N(ng)-1
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            buoy2(i,j,k)=bvf(i,j,k)
          END DO
        END DO
      END DO
# ifdef N2S2_HORAVG
!
!-----------------------------------------------------------------------
!  Smooth horizontally buoyancy and shear.  Use buoy2(:,:,0) and
!  shear2(:,:,0) as scratch utility array.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
            shear2(Istr-1,j,k)=shear2(Istr,j,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
            shear2(Iend+1,j,k)=shear2(Iend,j,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            shear2(i,Jstr-1,k)=shear2(i,Jstr,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            shear2(i,Jend+1,k)=shear2(i,Jend,k)
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          shear2(Istr-1,Jstr-1,k)=shear2(Istr,Jstr,k)
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          shear2(Istr-1,Jend+1,k)=shear2(Istr,Jend,k)
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          shear2(Iend+1,Jstr-1,k)=shear2(Iend,Jstr,k)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          shear2(Iend+1,Jend+1,k)=shear2(Iend,Jend,k)
        END IF
!
!  Average horizontally.
!
        DO j=Jstr-1,Jend
          DO i=Istr-1,Iend
            buoy2(i,j,0)=0.25_r8*(buoy2(i,j  ,k)+buoy2(i+1,j  ,k)+      &
     &                            buoy2(i,j+1,k)+buoy2(i+1,j+1,k))
            shear2(i,j,0)=0.25_r8*(shear2(i,j  ,k)+shear2(i+1,j  ,k)+   &
     &                             shear2(i,j+1,k)+shear2(i+1,j+1,k))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            buoy2(i,j,k)=0.25_r8*(buoy2(i,j  ,0)+buoy2(i-1,j  ,0)+      &
     &                            buoy2(i,j-1,0)+buoy2(i-1,j-1,0))
            shear2(i,j,k)=0.25_r8*(shear2(i,j  ,0)+shear2(i-1,j  ,0)+   &
     &                             shear2(i,j-1,0)+shear2(i-1,j-1,0))
          END DO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Time-step advective terms.
!-----------------------------------------------------------------------
!
!  At entry, it is assumed that the turbulent kinetic energy fields
!  "tke" and "gls", at time level "nnew", are set to its values at
!  time level "nstp" multiplied by the grid box thicknesses Hz
!  (from old time step and at W-points).
!
      DO k=1,N(ng)-1
# ifdef K_C2ADVECTION
!
!  Second-order, centered differences advection.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff=0.25_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            FXK(i,j)=cff*(tke(i,j,k,3)+tke(i-1,j,k,3))
            FXP(i,j)=cff*(gls(i,j,k,3)+gls(i-1,j,k,3))
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff=0.25_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            FEK(i,j)=cff*(tke(i,j,k,3)+tke(i,j-1,k,3))
            FEP(i,j)=cff*(gls(i,j,k,3)+gls(i,j-1,k,3))
          END DO
        END DO
# else
        DO j=Jstr,Jend
          DO i=Istrm1,Iendp2
            gradK(i,j)=(tke(i,j,k,3)-tke(i-1,j,k,3))
#  ifdef MASKING
            gradK(i,j)=gradK(i,j)*umask(i,j)
#  endif
            gradP(i,j)=(gls(i,j,k,3)-gls(i-1,j,k,3))
#  ifdef MASKING
            gradP(i,j)=gradP(i,j)*umask(i,j)
#  endif
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              gradK(Istr-1,j)=gradK(Istr,j)
              gradP(Istr-1,j)=gradP(Istr,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              gradK(Iend+2,j)=gradK(Iend+1,j)
              gradP(Iend+2,j)=gradP(Iend+1,j)
            END DO
          END IF
        END IF
#  ifdef K_C4ADVECTION
!
!  Fourth-order, centered differences advection.
!
        cff1=1.0_r8/6.0_r8
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff=0.5_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            FXK(i,j)=cff*0.5_r8*(tke(i-1,j,k,3)+tke(i,j,k,3)-           &
     &                           cff1*(gradK(i+1,j)-gradK(i-1,j)))
            FXP(i,j)=cff*0.5_r8*(gls(i-1,j,k,3)+gls(i,j,k,3)-           &
     &                           cff1*(gradP(i+1,j)-gradP(i-1,j)))
          END DO
        END DO
#  else
!
!  Third-order, upstream bias advection with velocity dependent
!  hyperdiffusion.
!
        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1
            curvK(i,j)=gradK(i+1,j)-gradK(i,j)
            curvP(i,j)=gradP(i+1,j)-gradP(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff=0.5_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            IF (cff.gt.0.0_r8) THEN
              cff1=curvK(i-1,j)
              cff2=curvP(i-1,j)
            ELSE
              cff1=curvK(i,j)
              cff2=curvP(i,j)
            END IF
            FXK(i,j)=cff*0.5_r8*(tke(i-1,j,k,3)+tke(i,j,k,3)-           &
     &                           Gadv*cff1)
            FXP(i,j)=cff*0.5_r8*(gls(i-1,j,k,3)+gls(i,j,k,3)-           &
     &                           Gadv*cff2)
          END DO
        END DO
#  endif
        DO j=Jstrm1,Jendp2
          DO i=Istr,Iend
            gradK(i,j)=(tke(i,j,k,3)-tke(i,j-1,k,3))
#  ifdef MASKING
            gradK(i,j)=gradK(i,j)*vmask(i,j)
#  endif
            gradP(i,j)=(gls(i,j,k,3)-gls(i,j-1,k,3))
#  ifdef MASKING
            gradP(i,j)=gradP(i,j)*vmask(i,j)
#  endif
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              gradK(i,Jstr-1)=gradK(i,Jstr)
              gradP(i,Jstr-1)=gradP(i,Jstr)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              gradK(i,Jend+2)=gradK(i,Jend+1)
              gradP(i,Jend+2)=gradP(i,Jend+1)
            END DO
          END IF
        END IF
#  ifdef K_C4ADVECTION
        cff1=1.0_r8/6.0_r8
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff=0.5_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            FEK(i,j)=cff*0.5_r8*(tke(i,j-1,k,3)+tke(i,j,k,3)-           &
     &                           cff1*(gradK(i,j+1)-gradK(i,j-1)))
            FEP(i,j)=cff*0.5_r8*(gls(i,j-1,k,3)+gls(i,j,k,3)-           &
     &                           cff1*(gradP(i,j+1)-gradP(i,j-1)))
          END DO
        END DO
#  else
        DO j=Jstr-1,Jend+1
          DO i=Istr,Iend
            curvK(i,j)=gradK(i,j+1)-gradK(i,j)
            curvP(i,j)=gradP(i,j+1)-gradP(i,j)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff=0.5_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            IF (cff.gt.0.0_r8) THEN
              cff1=curvK(i,j-1)
              cff2=curvP(i,j-1)
            ELSE
              cff1=curvK(i,j)
              cff2=curvP(i,j)
            END IF
            FEK(i,j)=cff*0.5_r8*(tke(i,j-1,k,3)+tke(i,j,k,3)-           &
     &                           Gadv*cff1)
            FEP(i,j)=cff*0.5_r8*(gls(i,j-1,k,3)+gls(i,j,k,3)-           &
     &                           Gadv*cff2)
          END DO
        END DO
#  endif
# endif
!
!  Time-step horizontal advection.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-                            &
     &                      cff*(FXK(i+1,j)-FXK(i,j)+                   &
     &                           FEK(i,j+1)-FEK(i,j))
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-                            &
     &                      cff*(FXP(i+1,j)-FXP(i,j)+                   &
     &                           FEP(i,j+1)-FEP(i,j))
          END DO
        END DO
      END DO
!
! Compute vertical advection.
!
      DO j=Jstr,Jend
# ifdef K_C2ADVECTION
        DO k=1,N(ng)
          DO i=Istr,Iend
            cff=0.25_r8*(W(i,j,k)+W(i,j,k-1))
            FCK(i,k)=cff*(tke(i,j,k,3)+tke(i,j,k-1,3))
            FCP(i,k)=cff*(gls(i,j,k,3)+gls(i,j,k-1,3))
          END DO
        END DO
# else
        cff1=7.0_r8/12.0_r8
        cff2=1.0_r8/12.0_r8
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            cff=0.5*(W(i,j,k)+W(i,j,k-1))
            FCK(i,k)=cff*(cff1*(tke(i,j,k-1,3)+                         &
     &                          tke(i,j,k  ,3))-                        &
     &                    cff2*(tke(i,j,k-2,3)+                         &
     &                          tke(i,j,k+1,3)))
            FCP(i,k)=cff*(cff1*(gls(i,j,k-1,3)+                         &
     &                          gls(i,j,k  ,3))-                        &
     &                    cff2*(gls(i,j,k-2,3)+                         &
     &                          gls(i,j,k+1,3)))
          END DO
        END DO
        cff1=1.0_r8/3.0_r8
        cff2=5.0_r8/6.0_r8
        cff3=1.0_r8/6.0_r8
         DO i=Istr,Iend
          cff=0.5_r8*(W(i,j,0)+W(i,j,1))
          FCK(i,1)=cff*(cff1*tke(i,j,0,3)+                              &
     &                  cff2*tke(i,j,1,3)-                              &
     &                  cff3*tke(i,j,2,3))
          FCP(i,1)=cff*(cff1*gls(i,j,0,3)+                              &
     &                  cff2*gls(i,j,1,3)-                              &
     &                  cff3*gls(i,j,2,3))
          cff=0.5_r8*(W(i,j,N(ng))+W(i,j,N(ng)-1))
          FCK(i,N(ng))=cff*(cff1*tke(i,j,N(ng)  ,3)+                    &
     &                      cff2*tke(i,j,N(ng)-1,3)-                    &
     &                      cff3*tke(i,j,N(ng)-2,3))
          FCP(i,N(ng))=cff*(cff1*gls(i,j,N(ng)  ,3)+                    &
     &                      cff2*gls(i,j,N(ng)-1,3)-                    &
     &                      cff3*gls(i,j,N(ng)-2,3))
        END DO
# endif
!
!  Time-step vertical advection term.
!
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-                            &
     &                      cff*(FCK(i,k+1)-FCK(i,k))
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-                            &
     &                      cff*(FCP(i,k+1)-FCP(i,k))
          END DO
        END DO
!
!----------------------------------------------------------------------
!  Compute vertical mixing, turbulent production and turbulent
!  dissipation terms.
!----------------------------------------------------------------------
!
!  Set term for vertical mixing of turbulent fields.
!
        cff=-0.5_r8*dt(ng)
        DO k=1,N(ng)
          DO i=Istr,Iend
            FCK(i,k)=cff*(Akk(i,j,k)+Akk(i,j,k-1))/Hz(i,j,k)
            CF(i,k)=0.0_r8
          END DO
        END DO
!
!  Compute production and dissipation terms.
!
        cff3=my_E2/(vonKar*vonKar)
        DO k=1,N(ng)-1
          DO i=Istr,Iend
!
!  Compute shear and bouyant production of turbulent energy (m3/s3)
!  at W-points (ignore small negative values of buoyancy).
!
            IF ((buoy2(i,j,k).gt.-5.0E-5_r8).and.                       &
     &          (buoy2(i,j,k).lt.0.0_r8)) THEN
              strat2=0.0_r8
            ELSE
              strat2=buoy2(i,j,k)
            END IF
            Qprod=shear2(i,j,k)*(Akv(i,j,k)-Akv_bak(ng))-               &
     &            strat2*(Akt(i,j,k,itemp)-Akt_bak(itemp,ng))
!
!  Recalculate old time-step unlimited length scale.
!
            Ls_unlmt=MAX(eps,                                           &
     &                   gls(i,j,k,nstp)/(MAX(tke(i,j,k,nstp),eps)))
!
!  Time-step production term.
!
            cff1=0.5_r8*(Hz(i,j,k)+Hz(i,j,k+1))
            tke(i,j,k,nnew)=tke(i,j,k,nnew)+                            &
     &                      dt(ng)*cff1*Qprod*2.0_r8
            gls(i,j,k,nnew)=gls(i,j,k,nnew)+                            &
     &                      dt(ng)*cff1*Qprod*my_E1*Ls_unlmt
!
!  Compute dissipation of turbulent energy (m3/s3).  Add in vertical
!  mixing term.
!
            Qdiss=dt(ng)*SQRT(tke(i,j,k,nstp))/(my_B1*Ls_unlmt)
            cff=Ls_unlmt*(1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,k))+           &
     &                    1.0_r8/(z_w(i,j,k)-z_w(i,j,0)))
            Wscale=1.0_r8+cff3*cff*cff
            BCK(i,k)=cff1*(1.0_r8+2.0_r8*Qdiss)-FCK(i,k)-FCK(i,k+1)
            BCP(i,k)=cff1*(1.0_r8+Wscale*Qdiss)-FCK(i,k)-FCK(i,k+1)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time-step dissipation and vertical diffusion terms implicitly.
!-----------------------------------------------------------------------
!
!  Set surface and bottom boundary conditions.
!
        DO i=Istr,Iend
          tke(i,j,N(ng),nnew)=my_B1p2o3*0.5_r8*                         &
     &                        SQRT((sustr(i,j)+sustr(i+1,j))**2+        &
     &                             (svstr(i,j)+svstr(i,j+1))**2)
          gls(i,j,N(ng),nnew)=0.0_r8
          tke(i,j,0,nnew)=my_B1p2o3*0.5_r8*                             &
     &                    SQRT((bustr(i,j)+bustr(i+1,j))**2+            &
     &                         (bvstr(i,j)+bvstr(i,j+1))**2)
          gls(i,j,0,nnew)=0.0_r8
        END DO
!
!  Solve tri-diagonal system for "tke".
!
        DO i=Istr,Iend
          cff=1.0_r8/BCK(i,N(ng)-1)
          CF(i,N(ng)-1)=cff*FCK(i,N(ng)-1)
          tke(i,j,N(ng)-1,nnew)=cff*(tke(i,j,N(ng)-1,nnew)-             &
     &                               FCK(i,N(ng))*tke(i,j,N(ng),nnew))
        END DO
        DO k=N(ng)-2,1,-1
          DO i=Istr,Iend
            cff=1.0_r8/(BCK(i,k)-CF(i,k+1)*FCK(i,k+1))
            CF(i,k)=cff*FCK(i,k)
            tke(i,j,k,nnew)=cff*(tke(i,j,k,nnew)-                       &
     &                           FCK(i,k+1)*tke(i,j,k+1,nnew))
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-CF(i,k)*tke(i,j,k-1,nnew)
          END DO
        END DO
!
!  Solve tri-diagonal system for "gls".
!
        DO i=Istr,Iend
          cff=1.0_r8/BCP(i,N(ng)-1)
          CF(i,N(ng)-1)=cff*FCK(i,N(ng)-1)
          gls(i,j,N(ng)-1,nnew)=cff*(gls(i,j,N(ng)-1,nnew)-             &
     &                               FCK(i,N(ng))*gls(i,j,N(ng),nnew))
        END DO
        DO k=N(ng)-2,1,-1
          DO i=Istr,Iend
            cff=1.0_r8/(BCP(i,k)-CF(i,k+1)*FCK(i,k+1))
            CF(i,k)=cff*FCK(i,k)
            gls(i,j,k,nnew)=cff*(gls(i,j,k,nnew)-                       &
     &                           FCK(i,k+1)*gls(i,j,k+1,nnew))
          END DO
        END DO
        DO k=1,N(ng)-1,+1
          DO i=Istr,Iend
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-CF(i,k)*gls(i,j,k-1,nnew)
          END DO
        END DO
!
!---------------------------------------------------------------------
!  Compute vertical mixing coefficients (m2/s).
!---------------------------------------------------------------------
!
        DO i=Istr,Iend
          qint=0.0_r8
          qzint=0.0_r8
          DO k=1,N(ng)-1
            if (tke(i,j,k,nnew).gt.0.0_r8) then
              cff=0.5_r8*(Hz(i,j,k)+Hz(i,j,k+1))
              q=sqrt(tke(i,j,k,nnew))
              qint=qint+q*cff
              qzint=qzint-q*(z_r(i,j,k)-z_w(i,j,N(ng)))*cff
            endif
          END DO
          if (tke(i,j,N(ng),nnew).gt.0.0_r8) then
            cff=Hz(i,j,k)
            q=sqrt(tke(i,j,N(ng),nnew))
            qint=qint+q*cff
            qzint=qzint-q*(z_r(i,j,k)-z_w(i,j,N(ng)))*cff
          endif
          lt(i)=0.23_r8*qzint/MAX(qint,1.0e-20)
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
!
!  Compute turbulent length scale (m).  The length scale is only
!  limited in the K-related calculations and not in QL production,
!  dissipation, wall-proximity, etc.
!
            tke(i,j,k,nnew)=MAX(tke(i,j,k,nnew),my_qmin)
            gls(i,j,k,nnew)=MAX(gls(i,j,k,nnew),my_qmin)
            if (tke(i,j,k,nnew).gt.0.0_r8) then
              q=sqrt(tke(i,j,k,nnew))
            else
              q=0.0_r8
            endif
            ls=-1.0_r8*vonKar*(z_w(i,j,k)-z_w(i,j,N(ng)))
            if (bvf(i,j,k).gt.0) then
              lbinv=sqrt(bvf(i,j,k))/MAX(q,1.0e-20)
            else
              lbinv=0.0_r8
            endif
            Ls_unlmt=1.0_r8/MAX(ls,1.0e-20)+1.0_r8/MAX(lt(i),1.0e-20)   &
     &               +lbinv
            dum=Ls_unlmt
            Ls_unlmt=1.0_r8/Ls_unlmt
            Ls_lmt=MIN(Ls_unlmt,                                        &
     &                 my_lmax*SQRT(tke(i,j,k,nnew)/                    &
     &                         (MAX(0.0_r8,buoy2(i,j,k))+eps)))
!
!  Compute Galperin et al. (1988) nondimensional stability function,
!  Gh.  Then, compute nondimensional stability functions for tracers
!  (Sh) and momentum (Sm).  The limit on length scale, sets the lower
!  limit on Gh.
!
            Gh=MIN(my_Gh0,-buoy2(i,j,k)*Ls_lmt*Ls_lmt/                    &
     &                    tke(i,j,k,nnew))
            dum=Ls_lmt*Ls_lmt/ tke(i,j,k,nnew)*shear2(i,j,k)
            Gm=MAX(Ls_lmt*Ls_lmt/tke(i,j,k,nnew)*shear2(i,j,k),0.0_r8)
            if (tke(i,j,k,nnew).gt.1.0d-8) then
               alpha_c=MIN(q/sqrt(tke(i,j,k,nnew)*2.0_r8),Ls_lmt)
            else
               alpha_c=MIN(cffm6,Ls_lmt)
            endif
            cff=alpha_c*alpha_c
            phi1=1.0_r8+cff*cffm1*Gh
            phi4=cff*cffm4*Gm
            phi5=1.0_r8+cff*cffm5*Gh
            phi2=phi5+cff*cffm2*Gh
            phi3=phi5+cff*cffm3*Gh
            cff=1.0_r8/(phi1*phi3+phi2*phi4)
            Sh=alpha_C*my_A2*(phi1+3.0_r8*my_C1*phi4)*cff
            Sm=alpha_C*my_A1*(phi2-3.0_r8*my_C1*phi3)*cff
            Sh=MAX(Sh,0.0_r8)
            Sm=MAX(Sm,0.0_r8)
            Sq=3.0*Sm
# if defined RR_MIXING || defined JK_MIXING
            Sh=MAX(Sh,1.0e-20_r8)
            SmSh(i,j,k)=Sm/Sh
# endif
!
!  Compute vertical mixing (m2/s) coefficients of momentum and
!  tracers.  Average ql over the two timesteps rather than using
!  the new Lscale and just averaging tke.
!
            ql=0.5_r8*(Ls_lmt*SQRT(tke(i,j,k,nnew))+                    &
     &                 Lscale(i,j,k)*SQRT(tke(i,j,k,nstp)))
            Akv(i,j,k)=ql*Sm
            Akv(i,j,k)=MAX(Akv_bak(ng),Akv(i,j,k))
            DO itrc=1,NAT
              Akt(i,j,k,itrc)=ql*Sh
              Akt(i,j,k,itrc)=MAX(Akt_bak(itrc,ng),Akt(i,j,k,itrc))
            END DO
!
!  Compute vertical mixing (m2/s) coefficient of turbulent kinetic
!  energy.  Use original formulation (Mellor and Yamada 1982;
!  Blumberg 1991; Kantha and Clayson 1994).
!
            Akk(i,j,k)=Akk_bak(ng)+ql*Sq
!
!  Save limited length scale.
!
            Lscale(i,j,k)=Ls_lmt
          END DO
        END DO
      END DO
!
!     Set the surface valueto the one below,
!               so it is not the background value
!
        DO i=Istr,Iend
          DO j=Jstr,Jend
            Akv(i,j,N(ng))=Akv(i,j,N(ng)-1);
            do itrc=1,NAT
              Akt(i,j,N(ng),itrc)=Akt(i,j,N(ng)-1,itrc);
            end do
        END DO
      END DO
!
!     Remove the spikes
!
      DO i=Istr,Iend
        DO j=Jstr,Jend
!
!       First determine the running average, then remove spikes
!                    ignore the ends
!         Akv first
          Do k=2,N(ng)-1
            Akdum(k)=0.25_r8*Akv(i,j,k-1)+0.5_r8*Akv(i,j,k)+     &
     &                            0.25_r8*Akv(i,j,k+1)

          end do
          Do k=2,N(ng)-1
            If (abs(Akv(i,j,k)-Akdum(k)).gt.0.5_r8*Akdum(k)) then
              Akv(i,j,k) = 0.5_r8*(Akv(i,j,k-1) + Akv(i,j,k+1))
            end if
          end do
!
!         now Akt
          Do itrc=1,NAT
            Do k=2,N(ng)-1
              Akdum(k) = 0.25_r8*Akt(i,j,k-1,itrc) +             &
     &                   0.5_r8*Akt(i,j,k,itrc) +                &
     &                   0.25_r8*Akt(i,j,k+1,itrc)

            end do
            Do k=2,N(ng)-1
              If (abs(Akt(i,j,k,itrc)-Akdum(k)).gt.              &
     &                          0.5_r8*Akdum(k)) then
                 Akt(i,j,k,itrc) = 0.5_r8*(Akt(i,j,k-1,itrc) +    &
     &                                     Akt(i,j,k+1,itrc))
              end if
            end do
          end do
        END DO
      END DO
# if defined RR_MIXING
!      write (6,*) (Akv(16,25,k),k=levbmix(ng),levsmix(ng),5)
!      write (6,*) 'lev',levsmix(ng),levbmix(ng),N(ng)
!--------------------------------------------------------------------
!
!   Replace mid-water column diffusivities with diffusivities based
!                        on shear and Ri   - Robertson
!                   or   shear cubed       - Shear3
!
!--------------------------------------------------------------------
      If ((levbmix(ng).lt.levsmix(ng)).and.(levbmix(ng).gt.1).and.     &
     &                   (levsmix(ng).lt.N(ng))) then
        cff= 7.5_r8*1.85E-6_r8
        DO j=Jstr,Jend
          DO i=Istr,Iend
            DO k=levbmix(ng),levsmix(ng)
#  if SHEAR3_MIXING
!              cff= 0.2_r8 * sqrt(abs(shear2(i,j,k)*100.0_r8))
              Akv(i,j,k)=Max(Akv_bak(ng),cff**3*SmSh(i,j,k),           &
     &                                   Akv(i,j,k))
              DO itrc=1,NAT
                Akt(i,j,k,itrc)=Max(Akt_bak(ng),cff**3*,               &
     &                                   Akt(i,j,k,itrc))
              END DO
#  else
              IF ((bvf(i,j,k).gt.-5.0E-5_r8).and.                       &
     &            (bvf(i,j,k).lt.0.0_r8)) THEN
                strat2=0.0_r8
              ELSE
                strat2=-bvf(i,j,k)
              ENDIF
              strat2=max(strat2,1.0e-4_r8)
              Riinv=max(shear2(i,j,k)*1._r8/strat2,0.0_r8)
           if ((i.eq.16).and.(j.eq.25)) then
            write (6,*) 'Riinv',Riinv,shear2(i,j,k),                    &
     &              shear2(i,j,k)*1._r8/strat2
           endif
              if (Riinv.gt.1000.) then
                 write(6,*)'Huge Riinv',i,j,k,Riinv,bvf(i,j,k),         &
      &                     shear2(i,j,k),strat2
              endif
!  Limit Inverse Richardson no to > molecular diff. and < 1 m^2/s
!              Riinv=max(Riinv,1.333333333333_r8)
!  Limit Inverse Richardson no to > Stillinger 1983 and < 1 m^2/s
              Riinv=max(Riinv,3.2_r8)
              Riinv=min(Riinv,7.2E5_r8)
              Akv(i,j,k)=Max(Akv_bak(ng),cff*Riinv*SmSh(i,j,k),         &
     &                                   Akv(i,j,k))
!              Akv(i,j,k)=Akv_bak(ng)+cff*Riinv
           if ((i.eq.16).and.(j.eq.25)) then
      Akdum(k)=(0.2_r8 * sqrt(abs(shear2(i,j,k)*1.0_r8)))**3
             endif
              DO itrc=1,NAT
                Akt(i,j,k,itrc)=Max(Akt_bak(itrc,ng),cff*Riinv,         &
     &                              Akt(i,j,k,itrc))
!                Akt(i,j,k,itrc)=Akt_bak(itrc,ng)+cff*Riinv
              END DO
#  endif
            END DO
          END DO
       END DO
      END IF
# elif defined JK_MIXING
!--------------------------------------------------------------------
!
!   Replace mid-water column diffusivities with diffusivities based
!                        on overturns   - following Kymak
!
!--------------------------------------------------------------------

      If ((levbmix(ng).lt.levsmix(ng)).and.(levbmix(ng).gt.1).and.      &
     &                   (levsmix(ng).lt.N(ng))) then
        DO j=Jstr,Jend
          DO i=Istr,Iend
            pdenmax=0.0_r8
!
!           Sort to determine the overturn distance
            Do k=1,N(ng)
              rho_dum(k) = pden(i,j,k)
              pdenmax=MAX(abs(pden(i,j,k)),pdenmax)
              rho_sort(k) = rho_dum(k)
            END DO

            IF (pdenmax.gt.0.00000001_r8) then
              call Bubble(rho_sort,N(ng))
              Do k=1,N(ng)
                rho_dums(k) = rho_sort(k)
              END DO
              Do k=levbmix(ng),levsmix(ng)
                kdum=1
                found = .false.
                Do while ((.not.found).and.(kdum.lt.N(ng)))
                  If (rho_dum(k).eq.rho_sort(kdum)) then
                    found=.true.
                    rho_sort(kdum)=0.0_r8
                  else
                    kdum=kdum+1
                  end if
                end do
                L_t(k) = z_r(i,j,kdum)-z_r(i,j,k)
              END DO
              DO k=levbmix(ng),levsmix(ng)
                diffdum = - (g / rho0) *(rho(i,j,k+1)-rho(i,j,k))/      &
     &                       (z_r(i,j,k+1)-z_r(i,j,k))
                If (diffdum.gt.0.0_r8) then
                  diffdum=MAX(diffdum,0.0000001_r8)
                  bvf_sort(k) = sqrt(diffdum)
                else
                  diffdum=-1.0_r8
                  bvf_sort(k)= -1.0_r8
                end if
              END DO
              DO k=levbmix(ng),levsmix(ng)
                if (bvf_sort(k).ge.0.0_r8) then
                  diffdum = 0.2_r8 * L_t(k) * L_t(k) * bvf_sort(k)
          if (L_t(k).gt.0.0_r8) then
          if ((i.eq.16).and.(j.eq.25).and.(mod(k,5).eq.1)) then
             write (6,*) 'dif L_t',i,j,k,diffdum, L_t(k),kdum,          &
      &            rho_dum(k),(rho_dums(kk),kk=kdum-2,kdum+2)
          end if
          end if
                  Akv(i,j,k)=Max(Akv(i,j,k),diffdum*SmSh(i,j,k),        &
     &                           Akv_bak(ng))
                  DO itrc=1,NAT
                    Akt(i,j,k,itrc)= Max(Akt(i,j,k,itrc),diffdum,       &
     &                                   Akt_bak(itrc,ng))
                  END DO
               end if
             END DO
           END IF
          END DO
        END DO
      END IF
      write (6,*) 'new',(Akv(16,25,k),k=levbmix(ng),levsmix(ng),5)
# endif
# if defined RR_MIXING || defined JK_MIXING || defined SHEAR3_MIXING
!
!     Smooth in the interior mixing scheme
      If ((levbmix(ng).lt.levsmix(ng)).and.(levbmix(ng).gt.1).and.     &
     &                   (levsmix(ng).lt.N(ng))) then
        DO j=Jstr,Jend
          DO i=Istr,Iend
            k=levbmix(ng)
            Akv(i,j,k)=0.5_r8*(Akv(i,j,k-1)+Akv(i,j,k+1))
            DO itrc=1,NAT
              Akt(i,j,k,itrc)= 0.5_r8*(Akt(i,j,k-1,itrc)+              &
     &                                 Akt(i,j,k+1,itrc))
            END DO
            k=levsmix(ng)
            Akv(i,j,k)=0.5_r8*(Akv(i,j,k-1)+Akv(i,j,k+1))
            DO itrc=1,NAT
              Akt(i,j,k,itrc)= 0.5_r8*(Akt(i,j,k-1,itrc)+              &
     &                                 Akt(i,j,k+1,itrc))
            END DO
          END DO
        END DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL tkebc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nnew, nstp,                                      &
     &                 gls, tke)
      CALL akbc_tile (ng, tile,                                         &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 Akv, Akt)
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          tke(:,:,:,nnew))
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          gls(:,:,:,nnew))
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          Akv)
        DO itrc=1,NAT
          CALL exchange_w3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 0, N(ng),         &
     &                            Akt(:,:,:,itrc))
        END DO
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tke(:,:,:,nnew),                              &
     &                    gls(:,:,:,nnew),                              &
     &                    Akv)
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng), 1, NAT,         &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akt)
# endif
      RETURN
      END SUBROUTINE NN_corstep_tile
!     Order subroutine
      Subroutine Order (p,q)
      real*8, intent(inout) :: p
      real*8, intent(inout) :: q
      real*8 :: temp
      if (p<q) then
        temp = p
        p = q
        q = temp
      end if
      return
      END SUBROUTINE ORDER

!     Bubble (A,B,nn)    largest to smallest
      Subroutine Bubble(A,nn)
      real*8, intent(inout) :: A(1:nn)
      integer, intent(in) :: nn
      integer :: i, j
      do i = 1, nn
        do j = nn, i+1, -1
          Call Order(A(j-1),A(j))
        end do
      end do
      return
      END SUBROUTINE BUBBLE
#endif
      END MODULE NN_corstep_mod
