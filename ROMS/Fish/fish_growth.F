#include "cppdefs.h"
      MODULE fish_growth_mod
#if defined NONLINEAR && defined NEMURO_SAN
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group         Mark Hadfield   !
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
     &                            egg_dur, egg_num,                     &
     &                            ysac_dur, ysac_num,                   &
     &                            larv_dur, larv_num,                   &
     &                            juv_dur, juv_num, suba_num,           &
     &                            fmortN, fmortS, fmortF, fmortP,       &
     &                            fish_count, fish_list, fishnodes)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
      USE mod_types
      USE mod_grid
      USE mod_biology
      USE mod_fish
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
      real(r8), intent(inout) :: egg_dur(Nfish(ng))
      real(r8), intent(inout) :: egg_num(Nfish(ng))
      real(r8), intent(inout) :: ysac_dur(Nfish(ng))
      real(r8), intent(inout) :: ysac_num(Nfish(ng))
      real(r8), intent(inout) :: larv_dur(Nfish(ng))
      real(r8), intent(inout) :: larv_num(Nfish(ng))
      real(r8), intent(inout) :: juv_dur(Nfish(ng))
      real(r8), intent(inout) :: juv_num(Nfish(ng))
      real(r8), intent(inout) :: suba_num(Nfish(ng))
      real(r8), intent(inout) :: fmortN(Nfish(ng))
      real(r8), intent(inout) :: fmortS(Nfish(ng))
      real(r8), intent(inout) :: fmortF(Nfish(ng))
      real(r8), intent(inout) :: fmortP(Nspecies(ng),Nfish(ng))
!
!  Local variable declarations.
!
      integer :: i1, i2, j1, j2, i, j, ii, jj
      integer :: ifish, ifid, isp

      type(fishnode), pointer :: thisfish

      real(r8) :: Fweight, Fworth, Fx, Fy, Fz
      real(r8) :: SZooC, LZooC, PZooC
      real(r8) :: FSZoo, FLZoo, FPZoo, FTotZoo
      real(r8) :: CSZoo, CLZoo, CPZoo, CTZoo
      real(r8) :: dtdays, Csmp, Egst, Excr, Resp, Actv, SDA
      real(r8) :: mmol2gww, gww2mmol, Cmax, Uswim
      real(r8) :: fac
      real(r8) :: Ndmort, Fdmort, Tdmort
!
!
!-----------------------------------------------------------------------
!  Compute feeding terms on zooplankton to update growth
!-----------------------------------------------------------------------
!
      dtdays=dt(ng)*sec2day
!
      mmol2gww=1.0_r8
      gww2mmol=1.0_r8/mmol2gww
!
! -----  PART 1: GROWTH   -----
!
      DO i=LBi+1,UBi-1
        DO j=LBj+1,UBj-1
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              isp = idfish(species(ifid))
              IF (MyThread(ifid).and.bounded(ifid).and.alive(ifid).and. &
     &            (lifestage(ifid).ge.if_juvenile)) THEN
                Fweight=bioenergy(ifwwt,ifid)
                Fworth=bioenergy(ifworth,ifid)
                Fx=track(ixgrd,itime,ifid)
                Fy=track(iygrd,itime,ifid)
! Swimming velocity based on fish weight (in cm/s)
                Uswim=a_A(isp,ng)*Fweight**b_A(isp,ng)

! JUVENILE, SUBADULT, OR ADULT FISH
! Max. Consumption (NEEDS TEMPERATURE DEPENDENCE)
                Cmax=a_C(isp,ng)*Fweight**(-b_C(isp,ng))
! Zoo Consumption
                SZooC=mmol2gww*track(NFV(ng)-NT(ng)+iSzoo,itime,ifid)
                LZooC=mmol2gww*track(NFV(ng)-NT(ng)+iLzoo,itime,ifid)
                PZooC=mmol2gww*track(NFV(ng)-NT(ng)+iPzoo,itime,ifid)
! To catch small negative zooplankton concentrations
                FSZoo=MAX(0.0_r8,SZooC)*ZSpref(isp,ng)/K_ZS(isp,ng)
                FLZoo=MAX(0.0_r8,LZooC)*ZLpref(isp,ng)/K_ZL(isp,ng)
                FPZoo=MAX(0.0_r8,PZooC)*ZPpref(isp,ng)/K_ZP(isp,ng)
                FTotZoo=FSZoo+FLZoo+FPZoo
! Total Consumption on Zooplankton
                CSZoo=Cmax*FSZoo*CAL_ZF(isp,ng)/(1.0_r8+FTotZoo)
                CLZoo=Cmax*FLZoo*CAL_ZF(isp,ng)/(1.0_r8+FTotZoo)
                CPZoo=Cmax*FPZoo*CAL_ZF(isp,ng)/(1.0_r8+FTotZoo)
                Csmp=CSZoo+CLZoo+CPZoo
! Egestion and Excretion
                Egst=a_F(isp,ng)*Csmp
                Excr=a_E(isp,ng)*(Csmp-Egst)
! Respiration
                Actv=exp(d_R(isp,ng)*Uswim)
                Resp=a_R(isp,ng)*Fweight**(-b_R(isp,ng))*Actv*5.258
! Specific Dynamic Action
                SDA=a_S(isp,ng)*(Csmp-Egst)
! Update Fish Weight (g wet wt)
                bioenergy(ifwwt,ifid)=Fweight*(1.0_r8+dtdays*           &
     &                                (Csmp-Egst-Excr-Resp-SDA))
!
! NEMURO sink terms for ZOO consumption and for DON egestion/excretion
! (will need to be converted from mmolN/day to mmolN/m3 in NEMURO)
! Recompute FFish
! Formulation for semi-implicit NEMURO
                fac=Cmax*Fweight*Fworth/(1.0_r8+FTotZoo)
                feedback(iSphy,ifid)=0.0_r8
                feedback(iLphy,ifid)=0.0_r8
                feedback(iSzoo,ifid)=fac*ZSpref(isp,ng)/K_ZS(isp,ng)
                feedback(iLzoo,ifid)=fac*ZLpref(isp,ng)/K_ZL(isp,ng)
                feedback(iPzoo,ifid)=fac*ZPpref(isp,ng)/K_ZP(isp,ng)
                feedback(iPON_,ifid)=0.0_r8
!                fac=Cmax*gww2mmol/(1.0_r8+FTotZoo+FTotFish)
!                feedback(iSzoo,ifid)=fac*FSZoo*Fweight*Fworth
!                feedback(iLzoo,ifid)=fac*FLZoo*Fweight*Fworth
!                feedback(iPzoo,ifid)=fac*FPZoo*Fweight*Fworth
!                feedback(iPON_,ifid)=gww2mmol*(Egst+Excr)*              &
!     &                               Fweight*Fworth/CAL_ZF(isp,ng)
!
! ADDITIONAL OUTPUT VARIABLES
! Fish P-value
                bioenergy(ifpval,ifid)=Csmp/Cmax
! Fish consumption
                bioenergy(ifcsmPS,ifid)=0.0_r8
                bioenergy(ifcsmPL,ifid)=0.0_r8
                bioenergy(ifcsmZS,ifid)=CSZoo
                bioenergy(ifcsmZL,ifid)=CLZoo
                bioenergy(ifcsmZP,ifid)=CPZoo
! Fish respiration
                bioenergy(ifresp,ifid)=Resp

              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
! -----   PART 2: MORTALITY   -----
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              IF (MyThread(ifid).and.bounded(ifid).and.alive(ifid).and. &
     &            (lifestage(ifid).ge.if_subadult)) THEN
                Fworth=bioenergy(ifworth,ifid)
                Ndmort=Nymort(isp,ng)/days_year      ! natural mortality rate per day
                Fdmort=Fymort(isp,ng)/days_year      ! fishing mortality rate per day
! Total Mortality (Natural + Fishing)
                Tdmort=Ndmort+Fdmort
! Update Fish Worth
                bioenergy(ifworth,ifid)=Fworth-                         &
     &                        (1.0_r8-exp(-Tdmort*dtdays))*Fworth
! Compute yield: yield=(Fm/(Fm+Nm))*(1-exp(-(Fm+Nm)))*Worth
                bioenergy(ifyield,ifid)=bioenergy(ifyield,ifid)+        &
     &            (Fdmort/Tdmort)*(1.0_r8-exp(-Tdmort*dtdays))*Fworth
! Natural and fishing mortality
                fmortN(ifid)=Ndmort
                fmortF(ifid)=Fdmort
! Starvation mortality (need to be computed somewhere else)
                fmortS(ifid)=0.0_r8
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
! -----   PART 3: EARLY LIFE HISTORY (EGG, YSAC, LARVA, JUV)  -----
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              IF (lifestage(ifid).le.if_juvenile) THEN
                Fweight=bioenergy(ifwwt,ifid)
                Fworth=bioenergy(ifworth,ifid)
! Juvenile (note: growth is done in PART 1 with subadult and adult) 
                IF (lifestage(ifid).eq.if_juvenile) THEN
                  Tdmort=0.03_r8
                  IF (iyesteryear(ng).ne.iyears(ng)) THEN
                    suba_num(ifid)=bioenergy(ifworth,ifid)
                  ELSE
                    juv_dur(ifid)=juv_dur(ifid)+dtdays
                  END IF
! Larva
                ELSE IF (lifestage(ifid).eq.if_larva) THEN
                  Fweight=Fweight+0.05_r8*dtdays
                  Tdmort=0.2_r8
                  IF (Fweight.ge.1.0_r8) THEN
                    lifestage(ifid)=if_juvenile
                    Fweight=0.5_r8*Fwwt0(isp,ng)
                    juv_dur(ifid)=0.0_r8
                    juv_num(ifid)=bioenergy(ifworth,ifid)
                  ELSE
                    larv_dur(ifid)=larv_dur(ifid)+dtdays
                  END IF
! Yolk Sac
                ELSE IF (lifestage(ifid).eq.if_yolksac) THEN
                  Fweight=Fweight+0.25_r8*dtdays
                  Tdmort=0.35_r8
                  IF (Fweight.ge.1.0_r8) THEN
                    lifestage(ifid)=if_larva
                    Fweight=0.0_r8
                    larv_dur(ifid)=0.0_r8
                    larv_num(ifid)=bioenergy(ifworth,ifid)
                  ELSE
                    ysac_dur(ifid)=ysac_dur(ifid)+dtdays
                  END IF
! Egg
                ELSE IF (lifestage(ifid).eq.if_egg) THEN
                  Fweight=Fweight+0.33_r8*dtdays
                  Tdmort=0.25_r8
                  IF (Fweight.ge.1.0_r8) THEN
                    lifestage(ifid)=if_yolksac
                    Fweight=0.0_r8
                    ysac_dur(ifid)=0.0_r8
                    ysac_num(ifid)=bioenergy(ifworth,ifid)
                  ELSE
                    egg_dur(ifid)=egg_dur(ifid)+dtdays
                  END IF
                END IF
                bioenergy(ifwwt,ifid)=Fweight
                bioenergy(ifworth,ifid)=Fworth-                         &
     &                        (1.0_r8-exp(-Tdmort*dtdays))*Fworth
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
! -----   PART 4: UPDATE FISH ATTRIBUTES   -----
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              IF (MyThread(ifid).and.bounded(ifid)) THEN
! From Kenny's code:
                IF (bioenergy(ifworth,ifid).lt.0.001_r8) THEN
                  alive(ifid)=.FALSE.
                END IF
                IF ((bioenergy(ifwwt,ifid).lt.0.001_r8).and.            &
     &              (lifestage(ifid).ge.if_juvenile)) THEN
                  alive(ifid)=.FALSE.
                END IF
                IF (alive(ifid)) THEN
! Update fish age and birthday
                  bioenergy(ifage,ifid)=bioenergy(ifage,ifid)+          &
     &                                        dtdays/days_year
                  IF (iyesteryear(ng).ne.iyears(ng)) THEN
                    bioenergy(ifbday,ifid)=bioenergy(ifbday,ifid)+1.0_r8
                  END IF
! Upgrade subadults to adult upon reaching maturity (Fage0)
                  IF ((lifestage(ifid).eq.if_subadult).and.             &
     &                (bioenergy(ifbday,ifid).ge.Fage0(isp,ng))) THEN
                    lifestage(ifid)=if_adult
                  END IF
                  IF ((lifestage(ifid).eq.if_juvenile).and.             &
     &                (bioenergy(ifbday,ifid).ge.1.0_r8)) THEN
                    lifestage(ifid)=if_subadult
                  END IF
                ELSE
                  bounded(ifid) = .FALSE.
                  bioenergy(ifwwt,ifid)=0.0_r8
                  bioenergy(ifworth,ifid)=0.0_r8
                  bioenergy(ifage,ifid)=0.0_r8
                  bioenergy(ifbday,ifid)=0.0_r8
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
