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
     &                            MyThread, bounded, track, feedback,   &
     &                            bioenergy, alive, eatfish, lifestage, &
     &                            fish_count, fish_list, fishnodes)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
      USE mod_types
      USE mod_grid
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj, itime
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)
      integer, intent(inout) :: lifestage(Nfish(ng))

      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Nfish(ng))

      logical, intent(in) :: MyThread(Nfish(ng))
      logical, intent(in) :: bounded(Nfish(ng))
      logical, intent(inout) :: alive(Nfish(ng))
      logical, intent(in) :: eatfish(Nfish(ng))

      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfish(ng))
      real(r8), intent(inout) :: feedback(NT(ng),Nfish(ng))
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfish(ng))
!
!  Local variable declarations.
!
      integer :: i1, i2, j1, j2, i, j
      integer :: ifish, ifid, ofish, ofid, iFwt

      type(fishnode), pointer :: thisfish
      type(fishnode), pointer :: otherfish

      real(r8) :: Fweight, Fworth, CellVol
      real(r8) :: oFweight, oFworth, FishVul, PopVul, BioVul
      real(r8) :: SZooC, LZooC, PZooC
      real(r8) :: FSZoo, FLZoo, FPZoo, FTZoo, FFish
      real(r8) :: CSZoo, CLZoo, CPZoo, CTZoo, CFish
      real(r8) :: dtdays, Csmp, Egst, Excr, Resp, Actv, SDA
      real(r8) :: mmol2g, g2mmol, CAL_ZF, Cmax, Uswim, fac
      real(r8) :: K_ZS, K_ZL, K_ZP, K_Fish, Fpref
      real(r8) :: a_C, b_C, a_F, a_E, a_R, b_R, d_R, a_A, b_A, a_S
      real(r8) :: Fpred, Nymort, Ndmort, Fymort, Fdmort, Pdmort
!
!-----------------------------------------------------------------------
!  Compute feeding terms on zooplankton to update growth
!-----------------------------------------------------------------------
!
      dtdays=dt(ng)*sec2day
!
      mmol2g=1.0_r8
      g2mmol=1.0_r8/mmol2g
! All parameters below should be species dependent
!      CAL_ZF=2580.0_r8/4460.0_r8           ! age 0-1
      CAL_ZF=2580.0_r8/7775.0_r8           ! age 2+ (averge)
      K_ZS=1.15_r8
      K_ZL=1.15_r8
      K_ZP=1.15_r8
      K_Fish=0.5_r8
      Fpref=1.0_r8
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
! Natural mortality should be age-dependent
      Nymort=0.3_r8               ! natural mortality rate per year
      Fymort=0.0_r8               ! fishing mortality rate per year
      Ndmort=Nymort/365.0_r8       ! natural mortality rate per day
      Fdmort=Fymort/365.0_r8       ! fishing mortality rate per day
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              IF (MyThread(ifid).and.bounded(ifid).and.alive(ifid)) THEN
                Fweight=bioenergy(ifwwt,ifid)
                Fworth=bioenergy(ifworth,ifid)
!
! ADULT FISH
                IF (lifestage(ifid).eq.if_adult) THEN
! Part 1. Fish growth
! Max. Consumption (NEEDS TEMPERATURE DEPENDENCE)
                  Cmax=a_C*Fweight**(-b_C)
! Zoo Consumption
                  SZooC=mmol2g*track(NFV(ng)-NT(ng)+iSzoo,itime,ifid)
                  LZooC=mmol2g*track(NFV(ng)-NT(ng)+iLzoo,itime,ifid)
                  PZooC=mmol2g*track(NFV(ng)-NT(ng)+iPzoo,itime,ifid)
                  FSZoo=SZooC*ZSpref(ng)/K_ZS
                  FLZoo=LZooC*ZLpref(ng)/K_ZL
                  FPZoo=PZooC*ZPpref(ng)/K_ZP
                  FTZoo=(FSZoo+FLZoo+FPZoo)*CAL_ZF
! Fish Consumption
! PROBLEM WITH 3D: Should the fish be allowed to see entire water column
! For now, fish_by_cell has contribution of entire water column
                  IF (eatfish(ifid)) THEN
                    otherfish => fish_list(i,j) % next
! FishVul, PopVul, and BioVul should be by species (e.g., BioVul(ispecies))
                    FishVul=0.0_r8
                    PopVul=0.0_r8
                    BioVul=0.0_r8
                    DO ofish=1,fish_count(i,j)
                      ofid = otherfish % fish
                      oFweight=bioenergy(ifwwt,ofid)
                      oFworth=bioenergy(ifworth,ofid)
                      IF (ifid.ne.ofid) THEN              ! No canibalism
                        IF (Fweight.gt.(2.0_r8*oFweight)) THEN
                          FishVul=FishVul+oFweight
                          PopVul=PopVul+oFworth
                          BioVul=BioVul+oFweight*oFworth
                        END IF
                      END IF 
                      otherfish => otherfish % next
                    END DO
! Should really divide by volume, but which Hz given that fish_by_cell is 2D?
                    CellVol=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
                    FishVul=FishVul*CellVol
                    PopVul=PopVul*CellVol
                    BioVul=BioVul*CellVol
                    FFish=BioVul*Fpref/K_Fish
                  ELSE
                    FFish=0.0_r8
                  END IF
                  CFish=Cmax*FFish/(1.0_r8+FTZoo+FFish)
! Must kill individual fish eaten to get correct consumption rates
! Approach: Distribute mortality evenly over vulnerable individuals in cell
                  fac=CFish*Fweight*Fworth*CellVol/BioVul
                  IF (fac.ge.1.0_r8) THEN
                    CFish=CFish*(1.0_r8-(fac-1.0_r8))
! Recompute FFish based on new CFish
                    FFish=(1.0_r8+FTZoo)/(Cmax/CFish-1.0_r8)
                    fac=0.99_r8
                  END IF 
                  IF (eatfish(ifid)) THEN
                    otherfish => fish_list(i,j) % next
                    DO ofish=1,fish_count(i,j)
                      ofid = otherfish % fish
                      oFweight=bioenergy(ifwwt,ofid)
                      IF (ifid.ne.ofid) THEN              ! No canibalism
                        IF (Fweight.gt.(2.0_r8*oFweight)) THEN
! bioenergy(ifeatme) needs to be reset to zero. In main3d or step_fish or here?
                          bioenergy(ifeatme,ofid)=                      &
     &                      bioenergy(ifeatme,ofid)+LOG(1.0_r8-fac)
                        END IF
                      END IF
                      otherfish => otherfish % next
                    END DO
                  END IF
! Total Consumption on Zooplankton and Fish
                  CSZoo=Cmax*FSZoo*CAL_ZF/(1.0_r8+FTZoo+FFish)
                  CLZoo=Cmax*FLZoo*CAL_ZF/(1.0_r8+FTZoo+FFish)
                  CPZoo=Cmax*FPZoo*CAL_ZF/(1.0_r8+FTZoo+FFish)
                  Csmp=CSZoo+CLZoo+CPZoo+CFish
! Fish P-value
                  bioenergy(ifpval,ifid)=Csmp/Cmax
! Egestion and Excretion
                  Egst=a_F*Csmp
                  Excr=a_E*(Csmp-Egst)
! (will need to be converted from mmolN/day to mmolN/m3 in NEMURO)
! Respiration (Uswim in cm/s)
                  Uswim=a_A*Fweight**b_A
                  Actv=exp(d_R*Uswim)
                  Resp=a_R*Fweight**(-b_R)*Actv*5.258
! Specific Dynamic Action
                  SDA=a_S*(Csmp-Egst)
! Update Fish Weight (g wet wt)
                  bioenergy(ifwwt,ifid)=Fweight*(1.0_r8+dtdays*         &
     &                                  (Csmp-Egst-Excr-Resp-SDA))
! From Kenny's code:
                  IF (bioenergy(ifwwt,ifid).lt.0.001_r8) THEN
                    bioenergy(ifwwt,ifid)=0.0_r8
                    alive(ifid)=.FALSE.
                  ENDIF
!
! Part 2. Fish Mortality
                  Fpred=1.0_r8                       ! No predatory fish
!                  Fpred=predbycell(i,j)/predavg     ! Compute from fish_by_cell
                  fac=1.0_r8-EXP(-Ndmort*0.1_r8)
                  fac=1.0_r8-fac*Fpred
                  fac=MIN(1.0_r8,MAX(1.0e-6_r8,fac))
                  Pdmort=-LOG(fac)
! Update Fish Worth
                  bioenergy(ifworth,ifid)=Fworth*(1.0_r8-dtdays*        &
     &                                    (Pdmort+Ndmort*0.9_r8+Fdmort))
! From Kenny's code:
                  IF (bioenergy(ifworth,ifid).lt.0.0001_r8) THEN
                    bioenergy(ifworth,ifid)=0.0_r8
                    alive(ifid)=.FALSE.
                  ENDIF
!
! Part 3. NEMURO sink terms for ZOO consumption and for DON egestion/excretion
! (will need to be converted from mmolN/day to mmolN/m3 in NEMURO)
! Recompute FFish
                  fac=Cmax*g2mmol/(FSZoo+FLZoo+FPZoo+FFish/CAL_ZF)
                  feedback(iSzoo,ifid)=fac*FSZoo*Fweight*Fworth
                  feedback(iLzoo,ifid)=fac*FLZoo*Fweight*Fworth
                  feedback(iPzoo,ifid)=fac*FPZoo*Fweight*Fworth
                  feedback(iPON_,ifid)=g2mmol*(Egst+Excr)*Fweight*Fworth
!
! EGG, YOLK SAC, AND FEEDING LARVA
                ELSE
! Larva
                  IF (lifestage(ifid).eq.if_larva) THEN
                    Fweight=Fweight+0.11_r8*dtdays
                    Fworth=Fworth-0.05_r8*dtdays*Fworth
                    IF (Fweight.ge.1.0_r8) THEN
                      lifestage(ifid)=if_adult
                      Fweight=Fwwt0(ng)
                    ENDIF
! Yolk Sac
                  ELSE IF (lifestage(ifid).eq.if_yolksac) THEN
                    Fweight=Fweight+0.1_r8*dtdays
                    Fworth=Fworth-0.05_r8*dtdays*Fworth
                    IF (Fweight.ge.1.0_r8) THEN
                      lifestage(ifid)=if_larva
                      Fweight=0.0_r8
                    ENDIF
! Egg
                  ELSE IF (lifestage(ifid).eq.if_egg) THEN
                    Fweight=Fweight+0.21_r8*dtdays
                    Fworth=Fworth-0.1_r8*dtdays*Fworth
                    IF (Fweight.ge.1.0_r8) THEN
                      lifestage(ifid)=if_yolksac
                      Fweight=0.0_r8
                    ENDIF
                  ENDIF
                  bioenergy(ifwwt,ifid)=Fweight
                  bioenergy(ifworth,ifid)=Fworth
!
                END IF
!
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE fish_growth
#endif
      END MODULE fish_growth_mod
