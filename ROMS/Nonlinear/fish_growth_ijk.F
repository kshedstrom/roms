#include "cppdefs.h"
#if defined NONLINEAR && defined FLOATS && defined NEMURO_SAN
      SUBROUTINE fish_growth_ijk (ng, Fweight, Uswim,                   &
# ifdef MASKING
     &                            Bmask,                                &
# endif
     &                            B, Fgrowth)
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
      integer, intent(in) :: ng
      real(r8), intent(inout) :: Fgrowth
      real(r8), intent(in) :: Uswim, Fweight
# ifdef MASKING
      real(r8), intent(in) :: Bmask
# endif
      real(r8), intent(in) :: B(1:NT(ng))
!
!  Local variable declarations.
!
      real(r8) :: SZooC, LZooC, PZooC, TZooC
      real(r8) :: FSZoo, FLZoo, FPZoo, FTZoo
      real(r8) :: dtdays, Csmp, Egst, Excr, Resp, Actv, SDA
      real(r8) :: mmol2g, g2mmol, CAL_ZF, C_max
      real(r8) :: K_ZS, K_ZL, K_ZP
      real(r8) :: a_C, b_C, a_F, a_E, a_R, b_R, d_R, a_S
!
!
!-----------------------------------------------------------------------
!  Compute swimming velocity based on food concentration
!-----------------------------------------------------------------------
!
      mmol2g=1.0_r8
      g2mmol=1.0_r8/mmol2g
!      CAL_ZF=2580.0_r8/4460.0_r8           ! age 0-1
      CAL_ZF=2580.0_r8/7775.0_r8           ! age 2+ (averge)
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
      a_S=0.175_r8
!
! COMSUMPTION
      C_max=a_C*Fweight**(-b_C)
# ifdef MASKING
      SZooC=mmol2g*Bmask*B(iSzoo)
      LZooC=mmol2g*Bmask*B(iLzoo)
      PZooC=mmol2g*Bmask*B(iPzoo)
# else
      SZooC=mmol2g*B(iSzoo)
      LZooC=mmol2g*B(iLzoo)
      PZooC=mmol2g*B(iPzoo)
# endif
      FSZoo=SZooC*ZSpref(ng)/K_ZS
      FLZoo=LZooC*ZLpref(ng)/K_ZL
      FPZoo=PZooC*ZPpref(ng)/K_ZP
      FTZoo=1.0_r8/(1.0_r8+FSZoo+FLZoo+FPZoo)
      Csmp=C_max*(FSZoo+FLZoo+FPZoo)*FTZoo*CAL_ZF
!
! EGESTION AND EXCRETION
      Egst=a_F*Csmp
      Excr=a_E*(Csmp-Egst)
!
! RESPIRATION (Uswim in cm/s)
      Actv=exp(d_R*Uswim)
      Resp=a_R*Fweight**(-b_R)*Actv*5.258
!
! SPECIFIC DYNAMIC ACTION
      SDA=a_S*(Csmp-Egst)
!
! GROWTH POTENTIAL at Ir,Jr,Kr
      Fgrowth=Csmp-Egst-Excr-Resp-SDA
!
#else
      SUBROUTINE fish_growth_ijk
#endif
      RETURN
      END SUBROUTINE fish_growth_ijk
