#include "cppdefs.h"
      MODULE fish_growth_mod
#if defined NONLINEAR && defined NEMURO_SAN
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

      implicit none

      PRIVATE
      PUBLIC  :: fish_growth

      CONTAINS
!
!***********************************************************************
      SUBROUTINE fish_growth (ng, LBi, UBi, LBj, UBj, itime,            &
     &                            MyThread, bounded, track,             &
     &                            feedback, bioenergy, alive,           &
     &                            fish_count, fish_list, fishnodes)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
      USE mod_types
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj, itime
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)

      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Nfish(ng))

      logical, intent(in) :: MyThread(Nfish(ng))
      logical, intent(in) :: bounded(Nfish(ng))
      logical, intent(in) :: alive(Nfish(ng))

      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfish(ng))
      real(r8), intent(inout) :: feedback(NT(ng),Nfish(ng))
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfish(ng))
!
!  Local variable declarations.
!
      integer :: Ir, Jr, Kr
      integer :: i1, i2, j1, j2, i, j, ifish, fid, iFwt

      type(fishnode), pointer :: thisfish

      real(r8) :: SZooC, LZooC, PZooC
      real(r8) :: FSZoo, FLZoo, FPZoo, FTZoo
      real(r8) :: dtdays, Csmp, Egst, Excr, Resp, Actv, SDA, Fweight
      real(r8) :: mmol2g, g2mmol, CAL_ZF, C_max, Uswim, fac
      real(r8) :: K_ZS, K_ZL, K_ZP
      real(r8) :: a_C, b_C, a_F, a_E, a_R, b_R, d_R, a_A, b_A, a_S
!
!-----------------------------------------------------------------------
!  Compute feeding terms on zooplankton to update growth
!-----------------------------------------------------------------------
!
      dtdays=dt(ng)*sec2day
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
      a_A=15.0_r8                ! a_A=3.9_r8
      b_A=0.13_r8
      a_S=0.175_r8
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          IF (fish_count(i,j).gt.0) THEN
            DO ifish=1,fish_count(i,j)
!
              thisfish => fish_list(i,j) % next
              fid = thisfish % fish
!
              IF (MyThread(fid).and.bounded(fid).and.alive(fid)) THEN
                Fweight=bioenergy(ifwwt,fid)
! OPTION 1: ADULT FISH
                IF (Fweight.ge.0.0_r8) THEN       ! Adult fish
! Zoo Consumption
                  C_max=a_C*Fweight**(-b_C)
                  SZooC=mmol2g*track(NFV(ng)-NT(ng)+iSzoo,itime,fid)
                  LZooC=mmol2g*track(NFV(ng)-NT(ng)+iLzoo,itime,fid)
                  PZooC=mmol2g*track(NFV(ng)-NT(ng)+iPzoo,itime,fid)
                  FSZoo=SZooC*ZSpref(ng)/K_ZS
                  FLZoo=LZooC*ZLpref(ng)/K_ZL
                  FPZoo=PZooC*ZPpref(ng)/K_ZP
                  FTZoo=1.0_r8/(1.0_r8+FSZoo+FLZoo+FPZoo)
                  Csmp=C_max*(FSZoo+FLZoo+FPZoo)*FTZoo*CAL_ZF
! Fish p-value
!            bioenergy(ifpvl,l)=(FSZoo+FLZoo+FPZoo)*FTZoo
! Fish predation sink terms for zooplankton in NEMURO
! (will need to be converted from mmolN/day to mmolN/m3 in NEMURO)
                  fac=Nfish(ng)*C_max*FTZoo*g2mmol
                  feedback(iSzoo,fid)=fac*FSZoo*Fweight
                  feedback(iLzoo,fid)=fac*FLZoo*Fweight
                  feedback(iPzoo,fid)=fac*FPZoo*Fweight
! Fish Consumption
! fish-eat-fish routine from KR
! Egestion and Excretion
                  Egst=a_F*Csmp
                  Excr=a_E*(Csmp-Egst)
! Fish egestion and excretion source terms for PON in NEMURO
! (will need to be converted from mmolN/day to mmolN/m3 in NEMURO)
                  feedback(iPON_,fid)=g2mmol*(Egst+Excr)*Fweight
! Respiration (Uswim in cm/s)
                  Uswim=a_A*Fweight**b_A
                  Actv=exp(d_R*Uswim)
                  Resp=a_R*Fweight**(-b_R)*Actv*5.258
! Specific Dynamic Action
                  SDA=a_S*(Csmp-Egst)
! Fish Growth (g wet wt)
                  bioenergy(ifwwt,fid)=Fweight*(1.0_r8+dtdays*          &
     &                                 (Csmp-Egst-Excr-Resp-SDA))
! OPTION 2: EGG, YOLK SAC, AND FEEDING LARVA
                ELSE
                  iFwt=-INT(Fweight)
                  IF ((iFwt.ge.1).AND.(iFwt.lt.2)) THEN
                    Fweight=ABS(Fweight)+0.21_r8*dtdays          ! Egg
                  ELSE IF ((iFwt.ge.2).AND.(iFwt.lt.3)) THEN
                    Fweight=ABS(Fweight)+0.1_r8*dtdays           ! Yolk sac
                  ELSE IF ((iFwt.ge.3).AND.(iFwt.lt.4)) THEN
                    Fweight=ABS(Fweight)+0.11_r8*dtdays          ! Larva
                  END IF
                  IF (Fweight.ge.4.0_r8) THEN
                    bioenergy(ifwwt,fid)=Fwwt0(ng)
                  ELSE
                    bioenergy(ifwwt,fid)=-Fweight
                  END IF
                END IF
              END IF
!
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE fish_growth
#endif
      END MODULE fish_growth_mod
