#include "cppdefs.h"
#if defined NONLINEAR && defined NEMURO_SAN
      SUBROUTINE fish_growth_ijk (ng, Ipos, Jpos, Kpos, nnew, isp, ils, &
     &                         Fweight, Flength, Ftemp, Uswim, Fgrowth)
!
!svn $Id$
!=================================================== Jerome Fiechter ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine interpolates requested field at the float trajectory   !
!  locations.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     Lstr       Starting float index to process.                      !
!     Lend       Ending   float index to process.                      !
!     itime      Floats time level to process.                         !
!     ifield     ID of field to compute.                               !
!     gtype      Grid type. If negative, interpolate floats slopes.    !
!     maskit     Should the field be masked? Ignored if Land/Sea       !
!                 masking is not active.                               !
!     nudg       Vertical random walk term to be added to the field.   !
!     pm         Inverse grid spacing (1/m) in the XI-direction.       !
!     pn         Inverse grid spacing (1/m) in the ETA-direction.      !
!     Hz         Vertical thicknesses (m).                             !
!     Amask      Field Land/Sea mask.                                  !
!     A          Field to interpolate from.                            !
!     fishthread   Float parallel thread bounded switch.                 !
!     bounded    Float grid bounded status switch.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     track      Interpolated field: track(ifield,itime,:).            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
      USE mod_ocean
      USE mod_grid
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Ipos, Jpos, Kpos, nnew, isp, ils
      real(r8), intent(in) :: Uswim, Fweight, Flength, Ftemp
      real(r8), intent(inout) :: Fgrowth
!
!  Local variable declarations.
!
      real(r8) :: Bsz, Blz, Bpz, Csz, Clz, Cpz, Csum
      real(r8) :: dtdays, Csmp, Resp, Actv
      real(r8) :: mmol2gww, gww2mmol, Cmax
      real(r8) :: tt5, t5, t4, tt7, t7, t6, gcta, gctb
      real(r8) :: fofT4C, fofT4R, pvalue, Aeff
!
!
!-----------------------------------------------------------------------
!  Compute growth conditions based on food concentration
!-----------------------------------------------------------------------
!
      mmol2gww=1.0_r8
      gww2mmol=1.0_r8/mmol2gww
!
! Max. Consumption
      tt5=1.0_r8/(te2(ils,isp,ng)-te1(ils,isp,ng))
      t5=tt5*log((xk2(ils,isp,ng)*(1.0_r8-xk1(ils,isp,ng)))/            &
     &           (xk4(ils,isp,ng)*(1.0_r8-xk3(ils,isp,ng))))
      t4=exp(t5*(Ftemp-te1(ils,isp,ng)))
      tt7=1.0_r8/(te4(ils,isp,ng)-te3(ils,isp,ng))
      t7=tt7*log((xk3(ils,isp,ng)*(1.0_r8-xk4(ils,isp,ng)))/            &
     &           (xk4(ils,isp,ng)*(1.0_r8-xk3(ils,isp,ng))))
      t6=exp(t7*(te4(ils,isp,ng)-Ftemp))
      gcta=xk1(ils,isp,ng)*t4/(1.0_r8+xk1(ils,isp,ng)*                  &
     &                                    (t4-1.0_r8))
      gctb=xk4(ils,isp,ng)*t6/(1.0_r8+xk4(ils,isp,ng)*                  &
     &                                    (t6-1.0_r8))
      fofT4C=gcta*gctb

      Cmax=fofT4C*a_C(ils,isp,ng)*Fweight**(-b_C(ils,isp,ng))

! Zoo concentrations
# ifdef MASKING
      Bsz=MAX(0.0_r8,OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iSzoo)*            &
     &                            GRID(ng)%rmask(Ipos,Jpos))
      Blz=MAX(0.0_r8,OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iLzoo)*            &
     &                            GRID(ng)%rmask(Ipos,Jpos))
      Bpz=MAX(0.0_r8,OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iPzoo)*            &
     &                            GRID(ng)%rmask(Ipos,Jpos))
# else
      Bsz=MAX(0.0_r8,OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iSzoo))
      Blz=MAX(0.0_r8,OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iLzoo))
      Bpz=MAX(0.0_r8,OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iPzoo))
# endif
! Prey vulnerability
      Csz=mmol2gww*Bsz*ZSpref(ils,isp,ng)/K_ZS(ils,isp,ng)
      Clz=mmol2gww*Blz*ZLpref(ils,isp,ng)/K_ZL(ils,isp,ng)
      Cpz=mmol2gww*Bpz*ZPpref(ils,isp,ng)/K_ZP(ils,isp,ng)
      Csum=Csz+Clz+Cpz
! Fish P-value
      pvalue=(Csz+Clz+Csz)/(1.0_r8+Csum)
      pvalue=MIN(pvalue,pvalmax(ils,isp,ng))
! Assimilation efficiency
      Aeff=MIN(a_AE(ils,isp,ng)*Fweight**b_AE(ils,isp,ng),              &
     &                                  AEmax(ils,isp,ng))
! Realized comsumption
      Csmp=Cmax*pvalue*Aeff
! Swimming velocity based on fish weight (in cm/s)
      IF (activity(ils,isp,ng).gt.0.0_r8) THEN
        Actv=activity(ils,isp,ng)
      ELSE
        Actv=exp(d_R(ils,isp,ng)*Fswim(ils,isp,ng)*                     &
     &                              0.1_r8*Flength)
        fofT4R=exp(cr(ils,isp,ng)*(Ftemp-tr(ils,isp,ng)))
        Resp=a_R(ils,isp,ng)*Fweight**(-b_R(ils,isp,ng))*               &
     &                               foft4R*Actv*5.285_r8
      END IF
      Fgrowth=Csmp-Resp
!
#else
      SUBROUTINE fish_growth_ijk
#endif
      RETURN
      END SUBROUTINE fish_growth_ijk
