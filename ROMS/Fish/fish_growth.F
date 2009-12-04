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
      SUBROUTINE fish_growth (ng, LBi, UBi, LBj, UBj, LBk, UBk, itime,  &
     &                            pm, pn,                               &
# ifdef SOLVE3D
     &                            Hz,                                   &
# endif
     &                            MyThread, bounded, track, feedback,   &
     &                            bioenergy, alive, eatfish,            &
     &                            species, lifestage,                   &
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
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj, LBk, UBk, itime
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)
      integer, intent(inout) :: species(Nfish(ng))
      integer, intent(inout) :: lifestage(Nfish(ng))

      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Nfish(ng))

      logical, intent(in) :: MyThread(Nfish(ng))
      logical, intent(inout) :: bounded(Nfish(ng))
      logical, intent(inout) :: alive(Nfish(ng))
      logical, intent(in) :: eatfish(Nfish(ng))

      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
# endif
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfish(ng))
      real(r8), intent(inout) :: feedback(NT(ng),Nfish(ng))
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfish(ng))
!
!  Local variable declarations.
!
      integer :: i1, i2, j1, j2, i, j
      integer :: ifish, ifid, ofish, ofid, iFwt, isp

      type(fishnode), pointer :: thisfish
      type(fishnode), pointer :: otherfish

      real(r8) :: Fweight, Fworth, Fx, Fy, Fz
      real(r8) :: oFweight, oFworth, oFx, oFy, oFz
      real(r8) :: FishVul, PopVul, BioVul, CellVol
      real(r8) :: SZooC, LZooC, PZooC
      real(r8) :: FSZoo, FLZoo, FPZoo, FTZoo, FFish
      real(r8) :: CSZoo, CLZoo, CPZoo, CTZoo, CFish
      real(r8) :: dtdays, Csmp, Egst, Excr, Resp, Actv, SDA
      real(r8) :: mmol2g, g2mmol, Cmax, Uswim, Fpred, fac, dist
      real(r8) :: Ndmort, Fdmort, Pdmort, Tdmort
!
!-----------------------------------------------------------------------
!  Compute feeding terms on zooplankton to update growth
!-----------------------------------------------------------------------
!
      dtdays=dt(ng)*sec2day
!
      mmol2g=1.0_r8
      g2mmol=1.0_r8/mmol2g
!
! Zero out eatme contributions from previous time step
      DO i=LBi,UBi
        DO j=LBj,UBj
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              IF (MyThread(ifid).and.bounded(ifid).and.alive(ifid)) THEN
                bioenergy(ifeatme,ifid)=0.0_r8
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
! -----  PART 1: GROWTH, INCLUDING FISH-EAT-FISH   -----
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              isp = idfish(species(ifid))
              IF (MyThread(ifid).and.bounded(ifid).and.alive(ifid)) THEN
                Fweight=bioenergy(ifwwt,ifid)
                Fworth=bioenergy(ifworth,ifid)
                Fx=track(ixgrd,itime,ifid)
                Fy=track(iygrd,itime,ifid)
! Swimming velocity based on fish weight (in cm/s)
                Uswim=a_A(isp,ng)*Fweight**b_A(isp,ng)

! ADULT FISH
                IF (lifestage(ifid).eq.if_adult) THEN
! Max. Consumption (NEEDS TEMPERATURE DEPENDENCE)
                  Cmax=a_C(isp,ng)*Fweight**(-b_C(isp,ng))
! Zoo Consumption
                  SZooC=mmol2g*track(NFV(ng)-NT(ng)+iSzoo,itime,ifid)
                  LZooC=mmol2g*track(NFV(ng)-NT(ng)+iLzoo,itime,ifid)
                  PZooC=mmol2g*track(NFV(ng)-NT(ng)+iPzoo,itime,ifid)
! To catch small negative zooplankton concentrations
                  FSZoo=MAX(0.0_r8,SZooC)*ZSpref(isp,ng)/K_ZS(isp,ng)
                  FLZoo=MAX(0.0_r8,LZooC)*ZLpref(isp,ng)/K_ZL(isp,ng)
                  FPZoo=MAX(0.0_r8,PZooC)*ZPpref(isp,ng)/K_ZP(isp,ng)
                  FTZoo=(FSZoo+FLZoo+FPZoo)*CAL_ZF(isp,ng)
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
                      IF (MyThread(ofid).and.bounded(ofid).and.         &
     &                                       alive(ofid)) THEN
                        oFweight=bioenergy(ifwwt,ofid)
                        oFworth=bioenergy(ifworth,ofid)
                        oFx=track(ixgrd,itime,ofid)
                        oFy=track(iygrd,itime,ofid)
                        dist=(((Fx-oFx)/pm(i,j))**2+                    &
     &                        ((Fy-oFy)/pn(i,j))**2)**0.5_r8
                        IF ((ifid.ne.ofid).and.                         &
                            (Fweight.gt.(2.0_r8*oFweight)).and.         &
     &                      (dist.lt.0.01_r8*Uswim*dt(ng))) THEN
                          FishVul=FishVul+oFweight
                          PopVul=PopVul+oFworth
                          BioVul=BioVul+oFweight*oFworth
                        END IF 
                      END IF 
                      otherfish => otherfish % next
                    END DO
! Should really divide by volume, but which Hz given that fish_by_cell is 2D?
                    CellVol=pm(i,j)*pn(i,j)
                    FishVul=FishVul*CellVol
                    PopVul=PopVul*CellVol
                    BioVul=BioVul*CellVol
                    FFish=BioVul*Fpref(isp,ng)/K_Fish(isp,ng)
                    CFish=Cmax*FFish/(1.0_r8+FTZoo+FFish)
! Must kill individual fish eaten to get correct consumption rates
! Approach: Distribute mortality evenly over vulnerable individuals in cell
! Kenny's Code. Not sure on the "(1-(fac-1)" reduction
!                    fac=CFish*Fweight*Fworth*CellVol/BioVul
!                    IF (fac.ge.1.0_r8) THEN
!                      CFish=CFish*(1.0_r8-(fac-1.0_r8))
! Recompute FFish based on new CFish
!                      FFish=(1.0_r8+FTZoo)/(Cmax/CFish-1.0_r8)
!                      fac=0.99_r8
!                    END IF 
! JF interpratation of Kenny's code. Could be totally wrong.
                    IF ((CFish*Fweight*Fworth).ge.(BioVul/CellVol)) THEN
                      CFish=CFish*(BioVul/CellVol)/(Fweight*Fworth)
                      FFish=(1.0_r8+FTZoo)/(Cmax/CFish-1.0_r8)
                      fac=0.99_r8
                    ELSE
                      fac=CFish*Fweight*Fworth*CellVol/BioVul
                    END IF
                    IF (eatfish(ifid)) THEN
                      otherfish => fish_list(i,j) % next
                      DO ofish=1,fish_count(i,j)
                        ofid = otherfish % fish
                        IF (MyThread(ofid).and.bounded(ofid).and.       &
     &                                         alive(ofid)) THEN
                          oFweight=bioenergy(ifwwt,ofid)
                          oFx=track(ixgrd,itime,ofid)
                          oFy=track(iygrd,itime,ofid)
                          dist=(((Fx-oFx)/pm(i,j))**2+                  &
     &                          ((Fy-oFy)/pn(i,j))**2)**0.5_r8
                          IF ((ifid.ne.ofid).and.                       &
                              (Fweight.gt.(2.0_r8*oFweight)).and.       &
     &                        (dist.lt.0.01_r8*Uswim*dt(ng))) THEN
                            bioenergy(ifeatme,ofid)=                    &
     &                        bioenergy(ifeatme,ofid)+LOG(1.0_r8-fac)
                          END IF
                        END IF
                        otherfish => otherfish % next
                      END DO
                    END IF
                  ELSE
                    FFish=0.0_r8
                    CFish=0.0_r8
                  END IF
! Total Consumption on Zooplankton and Fish
                  CSZoo=Cmax*FSZoo*CAL_ZF(isp,ng)/(1.0_r8+FTZoo+FFish)
                  CLZoo=Cmax*FLZoo*CAL_ZF(isp,ng)/(1.0_r8+FTZoo+FFish)
                  CPZoo=Cmax*FPZoo*CAL_ZF(isp,ng)/(1.0_r8+FTZoo+FFish)
                  Csmp=CSZoo+CLZoo+CPZoo+CFish
! Fish P-value
                  bioenergy(ifpval,ifid)=Csmp/Cmax
! Egestion and Excretion
                  Egst=a_F(isp,ng)*Csmp
                  Excr=a_E(isp,ng)*(Csmp-Egst)
! (will need to be converted from mmolN/day to mmolN/m3 in NEMURO)
! Respiration
                  Actv=exp(d_R(isp,ng)*Uswim)
                  Resp=a_R(isp,ng)*Fweight**(-b_R(isp,ng))*Actv*5.258
! Specific Dynamic Action
                  SDA=a_S(isp,ng)*(Csmp-Egst)
! Update Fish Weight (g wet wt)
                  bioenergy(ifwwt,ifid)=Fweight*(1.0_r8+dtdays*         &
     &                                  (Csmp-Egst-Excr-Resp-SDA))
! From Kenny's code:
                  IF (bioenergy(ifwwt,ifid).lt.0.001_r8) THEN
                    alive(ifid)=.FALSE.
                  ENDIF
!
! NEMURO sink terms for ZOO consumption and for DON egestion/excretion
! (will need to be converted from mmolN/day to mmolN/m3 in NEMURO)
! Recompute FFish
                  fac=Cmax*g2mmol/                                      &
     &                (FSZoo+FLZoo+FPZoo+FFish/CAL_ZF(isp,ng))
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
                      Fweight=Fwwt0(isp,ng)
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
                  bioenergy(ifpval,ifid)=0.0_r8
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
! -----   PART 2: MORTALITY, INCLUDING FISH-EAT-FISH   -----
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
                Ndmort=Nymort(isp,ng)/days_year      ! natural mortality rate per day
                Fdmort=Fymort(isp,ng)/days_year      ! fishing mortality rate per day
! Natural and Predation Mortality
                Fpred=1.0_r8                       ! No predatory fish
!                Fpred=predbycell(i,j)/predavg     ! Compute from fish_by_cell
                fac=1.0_r8-EXP(-Ndmort*0.1_r8)
                fac=1.0_r8-fac*Fpred
                fac=MIN(1.0_r8,MAX(1.0e-6_r8,fac))
                Pdmort=-LOG(fac)
! Total Mortality (Natural + Predation + Fishing + Fish-eat-Fish)
! (Pdmort, Ndmort, Fdmort are > 0, and eatme is < 0)
                Tdmort=Pdmort+Ndmort*0.9_r8+Fdmort-                     &
     &                 bioenergy(ifeatme,ifid)
! Update Fish Worth
                bioenergy(ifworth,ifid)=Fworth*(1.0_r8-dtdays*Tdmort)
! From Kenny's code:
                IF (bioenergy(ifworth,ifid).lt.0.0001_r8) THEN
                  alive(ifid)=.FALSE.
                ENDIF
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
! -----   PART 3: UPDATE FISH AGE   -----
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              IF (MyThread(ifid).and.bounded(ifid)) THEN
                IF (alive(ifid)) THEN
                  bioenergy(ifage,ifid)=bioenergy(ifage,ifid)+          &
     &                                  dtdays/days_year
                ELSE
                  bounded(ifid) = .FALSE.
                  bioenergy(ifwwt,ifid)=-99.0_r8
                  bioenergy(ifworth,ifid)=-99.0_r8
                  bioenergy(ifage,ifid)=-99.0_r8
                  lifestage(ifid)=0
                END IF
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
