#include "cppdefs.h"
#ifdef EOS_PERTURB
      MODULE rho_eos_perturb_mod
# ifdef SOLVE3D
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes  "in situ" density and other associated       !
!  quantitites as a function of potential temperature,  salinity,      !
!  and pressure from a polynomial expression (Jackett & McDougall,     !
!  1992). The polynomial expression was found from fitting to 248      !
!  values  in the  oceanographic  ranges of  salinity,  potential      !
!  temperature,  and pressure.  It  assumes no pressure variation      !
!  along geopotential surfaces, that is, depth (meters; negative)      !
!  and pressure (dbar; assumed negative here) are interchangeable.     !
!                                                                      !
!  Check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                        !
!                                                                      !
!     alpha = 2.1014611551470d-04 (1/Celsius)                          !
!     beta  = 7.2575037309946d-04 (1/PSU)                              !
!     gamma = 3.9684764511766d-06 (1/Pa)                               !
!     den   = 1050.3639165364     (kg/m3)                              !
!     den1  = 1028.2845117925     (kg/m3)                              !
!     sound = 1548.8815240223     (m/s)                                !
!     bulk  = 23786.056026320     (Pa)                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of     !
!    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.    !
!    and Oceanic Techn., vol. 12, pp. 381-389.                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: rho_eos_perturb
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rho_eos_perturb (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
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
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 14, __LINE__, __FILE__)
#  endif
      CALL rho_eos_perturb_tile (ng, tile,                              &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
#  ifdef VAR_RHO_2D
     &                   GRID(ng) % Hz,                                 &
#  endif
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   OCEAN(ng) % t,                                 &
#  ifdef VAR_RHO_2D
     &                   COUPLING(ng) % rhoA,                           &
     &                   COUPLING(ng) % rhoS,                           &
#  endif
#  ifdef BV_FREQUENCY
     &                   MIXING(ng) % bvf,                              &
#  endif
#  if defined LMD_SKPP    || defined LMD_BKPP         || \
     defined BULK_FLUXES || defined BALANCE_OPERATOR
     &                   MIXING(ng) % alpha,                            &
     &                   MIXING(ng) % beta,                             &
#   ifdef LMD_DDMIX
     &                   MIXING(ng) % alfaobeta,                        &
#   endif
#  endif
     &                   OCEAN(ng) % pden,                              &
     &                   OCEAN(ng) % rho)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 14, __LINE__, __FILE__)
#  endif

      RETURN
      END SUBROUTINE rho_eos_perturb

#  ifdef NONLIN_EOS
!
!***********************************************************************
      SUBROUTINE rho_eos_perturb_tile (ng, tile,                        &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
#   ifdef MASKING
     &                         rmask,                                   &
#   endif
#   ifdef VAR_RHO_2D
     &                         Hz,                                      &
#   endif
     &                         z_r, z_w, pm, pn, t,                     &
#   ifdef VAR_RHO_2D
     &                         rhoA, rhoS,                              &
#   endif
#   ifdef BV_FREQUENCY
     &                         bvf,                                     &
#   endif
#   if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
     &                         alpha, beta,                             &
#    ifdef LMD_DDMIX
     &                         alfaobeta,                               &
#    endif
#   endif
     &                         pden,                                    &
     &                         rho)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_stoperturb
      USE set_stoperturb
      USE exchange_2d_mod
      USE exchange_3d_mod
#   ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#   endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#   ifdef ASSUMED_SHAPE
#    ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#    endif
#    ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#    endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#    ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:,LBj:)
      real(r8), intent(out) :: rhoS(LBi:,LBj:)
#    endif
#    ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:,LBj:,0:)
#    endif
#    if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
      real(r8), intent(out) :: alpha(LBi:,LBj:)
      real(r8), intent(out) :: beta(LBi:,LBj:)
#     ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:,LBj:,0:)
#     endif
#    endif
      real(r8), intent(out) :: pden(LBi:,LBj:,:)
      real(r8), intent(out) :: rho(LBi:,LBj:,:)
#   else
#    ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#    endif
#    ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#    endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#    ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhoS(LBi:UBi,LBj:UBj)
#    endif
#    ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
#    endif
#    if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
      real(r8), intent(out) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: beta(LBi:UBi,LBj:UBj)
#     ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#     endif
#    endif
      real(r8), intent(out) :: pden(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: rho(LBi:UBi,LBj:UBj,N(ng))
#   endif
!
!  Local variable declarations.
!
      integer :: i, ised, itrc, j, k
      integer :: i1, i2, j1, j2, k1, k2

      real(r8) :: SedDen, Tp, Tpr10, Ts, Tt, sqrtTs
#   ifdef BV_FREQUENCY
      real(r8) :: bulk_dn, bulk_up, den_dn, den_up
#   endif
      real(r8) :: cff, cff1, cff2
      real(r8) :: deltaT, deltaTx, deltaTy, deltaTz
      !! *DEBUG
      real(r8) :: dTmax, dTmean, dSmax, dSmean
      !! DEBUG*
      real(r8) :: deltaS, deltaSx, deltaSy, deltaSz
      real(r8) :: DbulkDS_tmp,DbulkDT_tmp,Dden1DS_tmp,Dden1DT_tmp
      real(r8) :: bulk_tmp,bulk0_tmp,bulk1_tmp,bulk2_tmp
      real(r8) :: den_tmp,den1_tmp
      real(r8) :: den_pert1,den1_pert1,den_pert2,den1_pert2
      real(r8) :: clamp_dT, clamp_dS
#   ifdef EOS_TDERIVATIVE
      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Scof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Tcof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: wrk
#   endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk0
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den1

#   include "set_bounds.h"

!
!-----------------------------------------------------------------------
!  Get random perturbations
!-----------------------------------------------------------------------
!
      CALL get_pseudo2D_delta(ng, LBi, UBi, LBj, UBj,                   &
     &                   eos_pert_rcorr(ng),                            &
     &                   eos_pert_tcorr(ng),                            &
     &                   grid_mean_dx(ng),grid_mean_dy(ng),             &
     &                   PERTURB(ng) % deltax,                          &
     &                   var_perturb_dx(ng))

      CALL get_pseudo2D_delta(ng, LBi, UBi, LBj, UBj,                   &
     &                   eos_pert_rcorr(ng),                            &
     &                   eos_pert_tcorr(ng),                            &
     &                   grid_mean_dx(ng),grid_mean_dy(ng),             &
     &                   PERTURB(ng) % deltay,                          &
     &                   var_perturb_dy(ng))

      CALL get_pseudo2D_delta(ng, LBi, UBi, LBj, UBj,                   &
     &                   eos_pert_rcorr(ng),                            &
     &                   eos_pert_tcorr(ng),                            &
     &                   grid_mean_dx(ng),grid_mean_dy(ng),             &
     &                   PERTURB(ng) % deltaz,                          &
     &                   var_perturb_dz(ng))

      clamp_dT=0.6
      clamp_dS=0.2

      !! *DEBUG
      dTmax=0.0
      dSmax=0.0
      dTmean=0.0
      dSmean=0.0
      !! DEBUG*
      DO j=JstrT,JendT
        DO k=1,N(ng)
          DO i=IstrT,IendT
!
!-----------------------------------------------------------------------
!  Compute T and S gradients (centered)
!-----------------------------------------------------------------------
!
            IF (rmask(i,j) .EQ. 0) THEN
              deltaTx = 0
              deltaTy = 0
              deltaTz = 0
              deltaSx = 0
              deltaSy = 0
              deltaSz = 0
            ELSE
              ! If we are on the edges, not centered
              i1=MAX(i-1,IstrT)
              i2=MIN(i+1,IendT)
              ! If we are near the coast, not centered
              IF (rmask(i1,j) .EQ. 0) i1=i
              IF (rmask(i2,j) .EQ. 0) i2=i
              deltaTx = (t(i2,j,k,nrhs,itemp) - t(i1,j,k,nrhs,itemp))/  &
     &               MAX(i2-i1,1)
              deltaSx = (t(i2,j,k,nrhs,isalt) - t(i1,j,k,nrhs,isalt))/  &
     &               MAX(i2-i1,1)
              ! If we are on the edges, not centered
              j1=MAX(j-1,JstrT)
              j2=MIN(j+1,JendT)
              ! If we are near the coast, not centered
              IF (rmask(i,j1) .EQ. 0) j1=j
              IF (rmask(i,j2) .EQ. 0) j2=j
              deltaTy = (t(i,j2,k,nrhs,itemp) - t(i,j1,k,nrhs,itemp))/  &
     &               MAX(j2-j1,1)
              deltaSy = (t(i,j2,k,nrhs,isalt) - t(i,j1,k,nrhs,isalt))/  &
     &               MAX(j2-j1,1)
              ! If we are at surface or bottom, not centered
              k1=MAX(k-1,1)
              k2=MIN(k+1,N(ng))
              deltaTz = (t(i,j,k2,nrhs,itemp) - t(i,j,k1,nrhs,itemp))/  &
     &               (k2-k1)
              deltaSz = (t(i,j,k2,nrhs,isalt) - t(i,j,k1,nrhs,isalt))/  &
     &               (k2-k1)
            END IF
!
!-----------------------------------------------------------------------
!  Multiply by perturbation
!-----------------------------------------------------------------------
!
            ! Clamp dT and dS to avoid too high values
            IF (ABS(deltaTx) .GE. clamp_dT)                             &
     &              deltaTx=SIGN(clamp_dT,deltaTx)
            IF (ABS(deltaTy) .GE. clamp_dT)                             &
     &              deltaTy=SIGN(clamp_dT,deltaTy)
            IF (ABS(deltaTz) .GE. clamp_dT)                             &
     &              deltaTz=SIGN(clamp_dT,deltaTz)
            IF (ABS(deltaSx) .GE. clamp_dS)                             &
     &              deltaSx=SIGN(clamp_dS,deltaSx)
            IF (ABS(deltaSy) .GE. clamp_dS)                             &
     &              deltaSy=SIGN(clamp_dS,deltaSy)
            IF (ABS(deltaSz) .GE. clamp_dS)                             &
     &              deltaSz=SIGN(clamp_dS,deltaSz)
            deltaTx = deltaTx*(PERTURB(ng) % deltax(i,j))
            deltaTy = deltaTy*(PERTURB(ng) % deltay(i,j))
            deltaTz = deltaTz*(PERTURB(ng) % deltaz(i,j))
            deltaSx = deltaSx*(PERTURB(ng) % deltax(i,j))
            deltaSy = deltaSy*(PERTURB(ng) % deltay(i,j))
            deltaSz = deltaSz*(PERTURB(ng) % deltaz(i,j))
!
!-----------------------------------------------------------------------
!  Final stochastic perturbation
!-----------------------------------------------------------------------
!
            deltaT = deltaTx + deltaTy + deltaTz
            deltaS = deltaSx + deltaSy + deltaSz
            !! *DEBUG
            dTmax = MAX(ABS(deltaT),dTmax)
            dSmax = MAX(ABS(deltaS),dSmax)
            !! DEBUG*
!
!-----------------------------------------------------------------------
!  Compute density with perturbations on T and S
!-----------------------------------------------------------------------
!
            CALL nonlinear_rho(t(i,j,k,nrhs,itemp)+deltaT,              &
#   ifdef SALINITY
     &                         t(i,j,k,nrhs,isalt)+deltaS,              &
#   endif
     &                         z_r(i,j,k),                              &
#   ifdef SEDIMENT
     &                         t(i,j,k,nrhs,:),                         &
#   endif
#   ifdef EOS_TDERIVATIVE
     &                         Dden1DS_tmp,Dden1DT_tmp,                 &
     &                         DbulkDS_tmp,DbulkDT_tmp,                 &
#   endif
     &                         bulk_tmp,bulk0_tmp,                      &
     &                         bulk1_tmp,bulk2_tmp,                     &
     &                         den1_pert1,den_pert1)

            CALL nonlinear_rho(t(i,j,k,nrhs,itemp)-deltaT,              &
#   ifdef SALINITY
     &                         t(i,j,k,nrhs,isalt)-deltaS,              &
#   endif
     &                         z_r(i,j,k),                              &
#   ifdef SEDIMENT
     &                         t(i,j,k,nrhs,:),                         &
#   endif
#   ifdef EOS_TDERIVATIVE
     &                         Dden1DS_tmp,Dden1DT_tmp,                 &
     &                         DbulkDS_tmp,DbulkDT_tmp,                 &
#   endif
     &                         bulk_tmp,bulk0_tmp,                      &
     &                         bulk1_tmp,bulk2_tmp,                     &
     &                         den1_pert2,den_pert2)
!
!-----------------------------------------------------------------------
!  Compute Non linear density without perturbations
!-----------------------------------------------------------------------
!
            CALL nonlinear_rho(t(i,j,k,nrhs,itemp),                     &
#   ifdef SALINITY
     &                         t(i,j,k,nrhs,isalt),                     &
#   endif
     &                         z_r(i,j,k),                              &
#   ifdef SEDIMENT
     &                         t(i,j,k,nrhs,:),                         &
#   endif
#   ifdef EOS_TDERIVATIVE
     &                         Dden1DS(i,k),Dden1DT(i,k),               &
     &                         DbulkDS(i,k),DbulkDT(i,k),               &
#   endif
     &                         bulk(i,k),bulk0(i,k),                    &
     &                         bulk1(i,k),bulk2(i,k),                   &
     &                         den1_tmp,den_tmp)

            den(i,k)  = (den_pert1+den_pert2)/2
            den1(i,k) = (den1_pert1+den1_pert2)/2
!            den(i,k)  = den_tmp
!            den1(i,k) = den1_tmp

#   ifdef MASKING
            den(i,k)=den(i,k)*rmask(i,j)
#   endif
          END DO
        END DO

#   ifdef VAR_RHO_2D
!
!-----------------------------------------------------------------------
!  Compute vertical averaged density (rhoA) and density perturbation
!  (rhoS) used in barotropic pressure gradient.
!-----------------------------------------------------------------------
!
        DO i=IstrT,IendT
          cff1=den(i,N(ng))*Hz(i,j,N(ng))
          rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          rhoA(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrT,IendT
            cff1=den(i,k)*Hz(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)
            rhoA(i,j)=rhoA(i,j)+cff1
          END DO
        END DO
        cff2=1.0_r8/rho0
        DO i=IstrT,IendT
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          rhoA(i,j)=cff2*cff1*rhoA(i,j)
          rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j)
        END DO
#   endif

#   if defined BV_FREQUENCY
!
!-----------------------------------------------------------------------
!  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
!  and vertical W-points:
!
!                  bvf = - g/rho d(rho)/d(z).
!
!  The density anomaly difference is computed by lowering/rising the
!  water parcel above/below adiabatically at W-point depth "z_w".
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)-1
          DO i=IstrT,IendT
            bulk_up=bulk0(i,k+1)-                                       &
     &              z_w(i,j,k)*(bulk1(i,k+1)-                           &
     &                          bulk2(i,k+1)*z_w(i,j,k))
            bulk_dn=bulk0(i,k  )-                                       &
     &              z_w(i,j,k)*(bulk1(i,k  )-                           &
     &                          bulk2(i,k  )*z_w(i,j,k))
            cff1=1.0_r8/(bulk_up+0.1_r8*z_w(i,j,k))
            cff2=1.0_r8/(bulk_dn+0.1_r8*z_w(i,j,k))
            den_up=cff1*(den1(i,k+1)*bulk_up)
            den_dn=cff2*(den1(i,k  )*bulk_dn)
            bvf(i,j,k)=-g*(den_up-den_dn)/                              &
     &                 (0.5_r8*(den_up+den_dn)*                         &
     &                  (z_r(i,j,k+1)-z_r(i,j,k)))
          END DO
        END DO
        DO i=IstrT,IendT
!!        bvf(i,j,0)=bvf(i,j,1)
!!        bvf(i,j,N(ng))=bvf(i,j,N(ng)-1)
          bvf(i,j,0)=0.0_r8
          bvf(i,j,N(ng))=0.0_r8
        END DO
#   endif

#   if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
!
!-----------------------------------------------------------------------
!  Compute thermal expansion (1/Celsius) and saline contraction
!  (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
#    ifdef LMD_DDMIX
        DO k=1,N(ng)
#    else
        DO k=N(ng),N(ng)
#    endif
          DO i=IstrT,IendT
            Tpr10=0.1_r8*z_r(i,j,k)
!
!  Compute thermal expansion and saline contraction coefficients.
!
            cff=bulk(i,k)+Tpr10
            cff1=Tpr10*den1(i,k)
            cff2=bulk(i,k)*cff
            wrk(i,k)=(den(i,k)+1000.0_r8)*cff*cff
            Tcof(i,k)=-(DbulkDT(i,k)*cff1+                              &
     &                  Dden1DT(i,k)*cff2)
            Scof(i,k)= (DbulkDS(i,k)*cff1+                              &
     &                  Dden1DS(i,k)*cff2)
#    ifdef LMD_DDMIX
            alfaobeta(i,j,k)=Tcof(i,k)/Scof(i,k)
#    endif
          END DO
          IF (k.eq.N(ng)) THEN
            DO i=IstrT,IendT
              cff=1.0_r8/wrk(i,N(ng))
              alpha(i,j)=cff*Tcof(i,N(ng))
              beta (i,j)=cff*Scof(i,N(ng))
            END DO
          END IF
        END DO
#   endif
!
!-----------------------------------------------------------------------
!  Load "in situ" density anomaly (kg/m3 - 1000) and potential
!  density anomaly (kg/m3 - 1000) referenced to the surface into global
!  arrays. Notice that this is done in a separate (i,k) DO-loops to
!  facilitate the adjoint.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrT,IendT
            rho(i,j,k)=den(i,k)
            pden(i,j,k)=(den1(i,k)-1000.0_r8)
#   ifdef MASKING
            pden(i,j,k)=pden(i,j,k)*rmask(i,j)
#   endif
          END DO
        END DO
      END DO
      IF (Master) THEN
      WRITE(6,*) "RE:(rho_eos_pert) t,iic,dTmax,dSmax=",  &
     &         tile,iic(ng), dTmax, dSmax
      END IF
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rho)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          pden)

#   if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
#    ifdef LMD_DDMIX
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          alfaobeta)
#    endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          alpha)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          beta)
#   endif

#   ifdef VAR_RHO_2D
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoA)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoS)
#   endif

#   ifdef BV_FREQUENCY
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          bvf)
#   endif
      END IF

#   ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rho, pden)

#    if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
#     ifdef LMD_DDMIX
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    alfaobeta)
#     endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    alpha, beta)
#    endif

#    ifdef VAR_RHO_2D
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rhoA, rhoS)
#    endif

#    ifdef BV_FREQUENCY
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bvf)
#    endif
#   endif

      RETURN
      END SUBROUTINE rho_eos_perturb_tile






      SUBROUTINE nonlinear_rho(temp,                                    &
#   ifdef SALINITY
     &                         salt,                                    &
#   endif
     &                         zdepth,                                  &
#   ifdef SEDIMENT
     &                         Sedtracers,                              &
#   endif
#   ifdef EOS_TDERIVATIVE
     &                         Ddens1DS,Ddens1DT,                       &
     &                         DbulkDS,DbulkDT,                         &
#   endif
     &                         bulk,bulk0,bulk1,bulk2,                  &
     &                         dens1,dens)
!
!=======================================================================
!  Nonlinear equation of state.  Notice that this equation of state
!  is only valid for potential temperature range of -2C to 40C and
!  a salinity range of 0 PSU to 42 PSU.
!  Because of the pressure effect on freezing and the iceshelf, we've
!  extended temp limit to -2.5 (msd 4/25/06)
!=======================================================================
!
      USE mod_eoscoef
#   ifdef SEDIMENT
      USE mod_sediment
#   endif

!
!  Imported variable declarations.
!
      real(r8), intent(in)  :: temp     ! Local
#   ifdef SALINITY
      real(r8), intent(in)  :: salt     ! Local
#   endif
      real(r8), intent(in)  :: zdepth   ! Local
#   ifdef EOS_TDERIVATIVE
      real(r8), intent(out) :: Ddens1DS,Ddens1DT  ! Local d(dens1)/d(S,T)
      real(r8), intent(out) :: DbulkDS,DbulkDT    ! Local d(bulk)/d(S)
#   endif
      real(r8), intent(out) :: bulk,bulk0,bulk1,bulk2
      real(r8), intent(out) :: dens1,dens         ! Local density
!
!  Local variable declarations.
!
      real(r8) :: SedDen, Tp, Tpr10, Ts, Tt, sqrtTs
      real(r8) :: cff, cff1
      real(r8), dimension(0:9) :: C
#   ifdef EOS_TDERIVATIVE
      real(r8), dimension(0:9) :: dCdT(0:9)
#   endif
!
!  Check temperature and salinity lower values. Assign depth to the
!  pressure.
!
            Tt=MAX(-2.5_r8,temp)
            Tt=MIN(40.0_r8,Tt)
#   ifdef SALINITY
            Ts=MAX(0.0_r8,salt)
            Ts=MIN(100.0_r8,Ts)
            sqrtTs=SQRT(Ts)
#   else
            Ts=0.0_r8
            sqrtTs=0.0_r8
#   endif
            Tp=zdepth
            Tpr10=0.1_r8*Tp
!
!-----------------------------------------------------------------------
!  Compute density (kg/m3) at standard one atmosphere pressure.
!-----------------------------------------------------------------------
!
            C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
            C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
            C(2)=V00+Tt*(V01+Tt*V02)
#   ifdef EOS_TDERIVATIVE
!
            dCdT(0)=Q01+Tt*(2.0_r8*Q02+Tt*(3.0_r8*Q03+Tt*(4.0_r8*Q04+   &
     &                      Tt*5.0_r8*Q05)))
            dCdT(1)=U01+Tt*(2.0_r8*U02+Tt*(3.0_r8*U03+Tt*4.0_r8*U04))
            dCdT(2)=V01+Tt*2.0_r8*V02
#   endif
!
            dens1=C(0)+Ts*(C(1)+sqrtTs*C(2)+Ts*W00)

#   ifdef EOS_TDERIVATIVE
!
!  Compute d(dens1)/d(S) and d(dens1)/d(T) derivatives used in the
!  computation of thermal expansion and saline contraction
!  coefficients.
!
            Ddens1DS=C(1)+1.5_r8*C(2)*sqrtTs+2.0_r8*W00*Ts
            Ddens1DT=dCdT(0)+Ts*(dCdT(1)+sqrtTs*dCdT(2))
#   endif
!
!-----------------------------------------------------------------------
!  Compute secant bulk modulus.
!-----------------------------------------------------------------------
!
            C(3)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
            C(4)=B00+Tt*(B01+Tt*(B02+Tt*B03))
            C(5)=D00+Tt*(D01+Tt*D02)
            C(6)=E00+Tt*(E01+Tt*(E02+Tt*E03))
            C(7)=F00+Tt*(F01+Tt*F02)
            C(8)=G01+Tt*(G02+Tt*G03)
            C(9)=H00+Tt*(H01+Tt*H02)
#   ifdef EOS_TDERIVATIVE
!
            dCdT(3)=A01+Tt*(2.0_r8*A02+Tt*(3.0_r8*A03+Tt*4.0_r8*A04))
            dCdT(4)=B01+Tt*(2.0_r8*B02+Tt*3.0_r8*B03)
            dCdT(5)=D01+Tt*2.0_r8*D02
            dCdT(6)=E01+Tt*(2.0_r8*E02+Tt*3.0_r8*E03)
            dCdT(7)=F01+Tt*2.0_r8*F02
            dCdT(8)=G02+Tt*2.0_r8*G03
            dCdT(9)=H01+Tt*2.0_r8*H02
#   endif
!
            bulk0=C(3)+Ts*(C(4)+sqrtTs*C(5))
            bulk1=C(6)+Ts*(C(7)+sqrtTs*G00)
            bulk2=C(8)+Ts*C(9)
            bulk =bulk0-Tp*(bulk1-Tp*bulk2)

#   if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
!
!  Compute d(bulk)/d(S) and d(bulk)/d(T) derivatives used
!  in the computation of thermal expansion and saline contraction
!  coefficients.
!
            DbulkDS=C(4)+sqrtTs*1.5_r8*C(5)-                       &
     &                   Tp*(C(7)+sqrtTs*1.5_r8*G00-Tp*C(9))
            DbulkDT=dCdT(3)+Ts*(dCdT(4)+sqrtTs*dCdT(5))-           &
     &                   Tp*(dCdT(6)+Ts*dCdT(7)-                        &
     &                       Tp*(dCdT(8)+Ts*dCdT(9)))
#   endif
!
!-----------------------------------------------------------------------
!  Compute local "in situ" density anomaly (kg/m3 - 1000).
!-----------------------------------------------------------------------
!
            cff=1.0_r8/(bulk+Tpr10)
            dens=dens1*bulk*cff
#   if defined SEDIMENT && defined SED_DENS
            SedDen=0.0_r8
            DO ised=1,NST
              itrc=idsed(ised)
              cff1=1.0_r8/Srho(ised,ng)
              SedDen=SedDen+                                            &
     &               Sedtracers(itrc)*                                  &
     &               (Srho(ised,ng)-den)*cff1
            END DO
            dens=dens+SedDen
#   endif
            dens=dens-1000.0_r8

      END SUBROUTINE nonlinear_rho

#  endif
# endif

      END MODULE rho_eos_perturb_mod
#endif


