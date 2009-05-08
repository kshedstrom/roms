#include "cppdefs.h"
#if defined NONLINEAR && defined FLOATS && defined NEMURO_SAN
      SUBROUTINE fish_growth (ng, lflt, itime, track, feedback,         &
     &           bioenergy)
!
!svn $Id: interp_floats.F 294 2009-01-09 21:37:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group         Mark Hadfield   !
!    Licensed under a MIT/X style license             John M. Klinck   !
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
!     MyThread   Float parallel thread bounded switch.                 !
!     bounded    Float grid bounded status switch.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     track      Interpolated field: track(ifield,itime,:).            !
!     feedback   Feedback to ecosystem (NPZD) model.                   !
!     bioenergy  Bioenergetic fish fields.                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, lflt, itime

      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfloats(ng))
      real(r8), intent(inout) :: feedback(NT(ng),Nfloats(ng))
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfloats(ng))
!
!  Local variable declarations.
!
      integer :: Ir, Jr, Kr
      integer :: i1, i2, j1, j2, i, j

      real(r8) :: SZooC, LZooC, PZooC
      real(r8) :: FSZoo, FLZoo, FPZoo, FTZoo
      real(r8) :: dtdays, Csmp, Egst, Excr, Resp, Actv
      real(r8) :: mmol2g, g2mmol, CAL_ZF, C_max, Uswim, fac
      real(r8) :: K_ZS, K_ZL, K_ZP
      real(r8) :: a_C, b_C, a_F, a_E, a_R, b_R, d_R, a_A, b_A
!      real(r8) :: Nfish
!
!
!-----------------------------------------------------------------------
!  Compute feeding terms on zooplankton to update growth
!-----------------------------------------------------------------------
!
      dtdays=dt(ng)*sec2day
!
!      Nfish=1.0_r8
!
      mmol2g=1.0_r8
      g2mmol=1.0_r8/mmol2g
      CAL_ZF=3800.0_r8/4460.0_r8
      K_ZS=1.15_r8
      K_ZL=1.15_r8
      K_ZP=1.15_r8
      a_C=0.642_r8
      b_C=0.256_r8
      a_F=0.16_r8
      a_E=0.10_r8
      a_R=0.0033_r8
      b_R=0.227_r8
      d_R=0.03_r8
      a_A=3.9_r8
      b_A=0.13_r8
!
! COMSUMPTION
      C_max=a_C*bioenergy(ifwwt,lflt)**(-b_C)
      SZooC=mmol2g*track(NFV(ng)-NT(ng)+iSzoo,itime,lflt)
      LZooC=mmol2g*track(NFV(ng)-NT(ng)+iLzoo,itime,lflt)
      PZooC=mmol2g*track(NFV(ng)-NT(ng)+iPzoo,itime,lflt)
      FSZoo=SZooC*ZSpref(ng)/K_ZS
      FLZoo=LZooC*ZLpref(ng)/K_ZL
      FPZoo=PZooC*ZPpref(ng)/K_ZP
      FTZoo=1.0_r8/(1.0_r8+FSZoo+FLZoo+FPZoo)
      Csmp=C_max*(FSZoo+FLZoo+FPZoo)*FTZoo*CAL_ZF
! Fish predation sink terms on zooplankton (mmolN/day) for NEMURO
! (will need to be converted to mmolN/m3 in NEMURO)
      fac=Nfish(ng)*C_max*FTZoo*g2mmol
      feedback(iSzoo,lflt)=fac*FSZoo*bioenergy(ifwwt,lflt)
      feedback(iLzoo,lflt)=fac*FLZoo*bioenergy(ifwwt,lflt)
      feedback(iPzoo,lflt)=fac*FPZoo*bioenergy(ifwwt,lflt)
!
! EGESTION AND EXCRETION
      Egst=a_F*Csmp
      Excr=a_E*(Csmp-Egst)
!
! RESPIRATION (Uswim in cm/s)
      Uswim=a_A*bioenergy(ifwwt,lflt)**b_A
      Actv=exp(d_R*Uswim)
      Resp=a_R*bioenergy(ifwwt,lflt)**(-b_R)*Actv*5.258
! Fish growth (g wet wt)
      bioenergy(ifwwt,lflt)=bioenergy(ifwwt,lflt)+                      &
     &                        dtdays*(Csmp-Egst-Excr-Resp)*             &
     &                        bioenergy(ifwwt,lflt)
!
      feedback(NT(ng),lflt)=bioenergy(ifwwt,lflt)
!
#else
      SUBROUTINE fish_growth
#endif
      RETURN
      END SUBROUTINE fish_growth
