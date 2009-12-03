#include "cppdefs.h"
#if defined NONLINEAR && defined NEMURO_SAN
      SUBROUTINE fish_growth_ijk (ng, Ipos, Jpos, Kpos, nnew,           &
     &                                isp, Fweight, Uswim, Fgrowth)
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
      USE mod_ocean
      USE mod_grid
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Ipos, Jpos, Kpos, nnew, isp
      real(r8), intent(in) :: Uswim, Fweight
      real(r8), intent(inout) :: Fgrowth
!
!  Local variable declarations.
!
      real(r8) :: SZooC, LZooC, PZooC, TZooC
      real(r8) :: FSZoo, FLZoo, FPZoo, FTZoo
      real(r8) :: dtdays, Csmp, Egst, Excr, Resp, Actv, SDA
      real(r8) :: mmol2g, g2mmol, C_max
!      real(r8) :: K_ZS, K_ZL, K_ZP, CAL_ZF
!      real(r8) :: a_C, b_C, a_F, a_E, a_R, b_R, d_R, a_S
!
!
!-----------------------------------------------------------------------
!  Compute growth conditions based on food concentration
!-----------------------------------------------------------------------
!
      mmol2g=1.0_r8
      g2mmol=1.0_r8/mmol2g
!
! Comsumption
      C_max=a_C(isp,ng)*Fweight**(-b_C(isp,ng))
# ifdef MASKING
      SZooC=mmol2g*GRID(ng)%rmask(Ipos,Jpos)*                           &
     &             OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iSzoo)
      LZooC=mmol2g*GRID(ng)%rmask(Ipos,Jpos)*                           &
     &             OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iLzoo)
      PZooC=mmol2g*GRID(ng)%rmask(Ipos,Jpos)*                           &
     &             OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iPzoo)
# else
      SZooC=mmol2g*OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iSzoo)
      LZooC=mmol2g*OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iLzoo)
      PZooC=mmol2g*OCEAN(ng)%t(Ipos,Jpos,Kpos,nnew,iPzoo)
# endif
      FSZoo=SZooC*ZSpref(isp,ng)/K_ZS(isp,ng)
      FLZoo=LZooC*ZLpref(isp,ng)/K_ZL(isp,ng)
      FPZoo=PZooC*ZPpref(isp,ng)/K_ZP(isp,ng)
      FTZoo=1.0_r8/(1.0_r8+FSZoo+FLZoo+FPZoo)
      Csmp=C_max*(FSZoo+FLZoo+FPZoo)*FTZoo*CAL_ZF(isp,ng)
! Egestion and Excretion
      Egst=a_F(isp,ng)*Csmp
      Excr=a_E(isp,ng)*(Csmp-Egst)
! Respiration (Uswim in cm/s)
      Actv=exp(d_R(isp,ng)*Uswim)
      Resp=a_R(isp,ng)*Fweight**(-b_R(isp,ng))*Actv*5.258
! Specific Dynamic Action
      SDA=a_S(isp,ng)*(Csmp-Egst)
! Growth Potential at Ipos,Jpos,Kpos
      Fgrowth=Csmp-Egst-Excr-Resp-SDA
!
#else
      SUBROUTINE fish_growth_ijk
#endif
      RETURN
      END SUBROUTINE fish_growth_ijk
