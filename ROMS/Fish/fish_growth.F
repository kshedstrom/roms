#include "cppdefs.h"
      MODULE fish_growth_mod
#if defined NONLINEAR && defined NEMURO_SAN
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group         Mark Hadfield   !
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
!     fishthread   Float parallel thread bounded switch.               !
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
      SUBROUTINE fish_growth (ng, tile, LBi, UBi, LBj, UBj, LBk, UBk,   &
     &                            itime, nnew, pm, pn,                  &
# ifdef SOLVE3D
     &                            Hz,                                   &
# endif
     &                            fishthread, bounded, track, feedback, &
     &                            bioenergy, alive, species,            &
     &                            lifestage,                            &
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
      integer, intent(in) :: ng, tile, itime, nnew
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)
      integer, intent(inout) :: species(Nfish(ng))
      integer, intent(inout) :: lifestage(Nfish(ng))

      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Nfish(ng))

      logical, intent(in) :: fishthread(Nfish(ng))
      logical, intent(inout) :: bounded(Nfish(ng))
      logical, intent(inout) :: alive(Nfish(ng))

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
      integer :: i1, i2, j1, j2, i, j, ii, jj, k
      integer :: ifish, ifid, isp, ifsp, ils

      integer :: Fcount(N(ng))

      type(fishnode), pointer :: thisfish

      real(r8) :: Fweight, Fworth, Flength, Ftemp, Fbiomass
      real(r8) :: Fx, Fy, Fz, Vcell, Bsz, Blz, Bpz
      real(r8) :: dtdays, dtsecs, Csmp, Resp, Actv, Cal_ZF, Csum
      real(r8) :: mmol2gww, gww2mmol, fofT4C, fofT4R, fac
      real(r8) :: tt5, t5, t4, tt7, t7, t6, gcta, gctb
      real(r8) :: pvalue, Aeff, LfromW, WfromL, deltaL1, deltaL2
      real(r8) :: dy_time, t_ramp, tday, fishJ, preyJ

      real(r8) :: BszAvg(N(ng))
      real(r8) :: BlzAvg(N(ng))
      real(r8) :: BpzAvg(N(ng))
      real(r8) :: CszTot(N(ng))
      real(r8) :: ClzTot(N(ng))
      real(r8) :: CpzTot(N(ng))
      real(r8) :: Csz_fac(N(ng))
      real(r8) :: Clz_fac(N(ng))
      real(r8) :: Cpz_fac(N(ng))
      real(r8) :: Csz(Nfish(ng))
      real(r8) :: Clz(Nfish(ng))
      real(r8) :: Cpz(Nfish(ng))
      real(r8) :: CmaxF(Nfish(ng))

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute growth terms for fish (larva, juvenile, adult)
!-----------------------------------------------------------------------
!
      dtsecs=3600.0_r8    ! one hour
      dtdays=dtsecs*sec2day
! dy_time: time in days, modulo days_year
! Note: only works if time origin is at the start of a year
      dy_time=REAL(INT(time(ng)/86400.0_r8/days_year))
      dy_time=time(ng)/86400.0_r8-days_year*dy_time
!
! mmolN/m3*0.001mol/mmol*14gdwN/mol*1gdwprey/0.07gdwN*5gwwprey/gdwprey=gwwprey/m3
      mmol2gww=1.0_r8
      gww2mmol=1.0_r8/mmol2gww
!
!      DO i=LBi+2,UBi-2
!        DO j=LBj+2,UBj-2
      DO i=Istr,Iend
        DO j=Jstr,Jend
          DO k=1,N(ng)
            BszAvg(k)=0.0_r8
            BlzAvg(k)=0.0_r8
            BpzAvg(k)=0.0_r8
            Fcount(k)=0
            CszTot(k)=0.0_r8
            ClzTot(k)=0.0_r8
            CpzTot(k)=0.0_r8
            Csz_fac(k)=1.0_r8
            Clz_fac(k)=1.0_r8
            Cpz_fac(k)=1.0_r8
          END DO
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              isp = idfish(species(ifid))
              IF (fishthread(ifid).and.bounded(ifid).and.               &
     &            alive(ifid).and. &
     &            (lifestage(ifid).ge.if_larva)) THEN
                Fweight=MAX(0.0_r8,bioenergy(ifwwt,ifid))
                Fworth=bioenergy(ifworth,ifid)
                Flength=bioenergy(iflngth,ifid)
                Fbiomass=Fweight*Fworth
                ils=lifestage(ifid)
                k=MIN(MAX(NINT(track(izgrd,itime,ifid)),1),N(ng))
                Csz(ifid)=0.0_r8
                Clz(ifid)=0.0_r8
                Cpz(ifid)=0.0_r8
                CmaxF(ifid)=0.0_r8
! Max. Consumption
                Ftemp=track(itemp+NFV(ng)-NT(ng),itime,ifid)
                tt5=1.0_r8/(te2(ils,isp,ng)-te1(ils,isp,ng))
                t5=tt5*LOG((xk2(ils,isp,ng)*(1.0_r8-xk1(ils,isp,ng)))/  &
     &                     (xk4(ils,isp,ng)*(1.0_r8-xk3(ils,isp,ng))))
                t4=EXP(t5*(Ftemp-te1(ils,isp,ng)))
                tt7=1.0_r8/(te4(ils,isp,ng)-te3(ils,isp,ng))
                t7=tt7*LOG((xk3(ils,isp,ng)*(1.0_r8-xk4(ils,isp,ng)))/  &
     &                     (xk4(ils,isp,ng)*(1.0_r8-xk3(ils,isp,ng))))
                t6=EXP(t7*(te4(ils,isp,ng)-Ftemp))
                gcta=xk1(ils,isp,ng)*t4/(1.0_r8+xk1(ils,isp,ng)*        &
     &                                              (t4-1.0_r8))
                gctb=xk4(ils,isp,ng)*t6/(1.0_r8+xk4(ils,isp,ng)*        &
     &                                              (t6-1.0_r8))
                fofT4C=gcta*gctb

                CmaxF(ifid)=fofT4C*a_C(ils,isp,ng)*Fweight**            &
     &                                    (-b_C(ils,isp,ng))
! Zoo concentrations
                Bsz=MAX(0.0_r8,track(NFV(ng)-NT(ng)+iSzoo,itime,ifid))
                Blz=MAX(0.0_r8,track(NFV(ng)-NT(ng)+iLzoo,itime,ifid))
                Bpz=MAX(0.0_r8,track(NFV(ng)-NT(ng)+iPzoo,itime,ifid))
! Keep track of average concentration in cell for overgrazing
                BszAvg(k)=BszAvg(k)+Bsz
                BlzAvg(k)=BlzAvg(k)+Blz
                BpzAvg(k)=BpzAvg(k)+Bpz
                Fcount(k)=Fcount(k)+1
! Prey vulnerability
                Csz(ifid)=mmol2gww*Bsz*ZSpref(ils,isp,ng)/              &
     &                                    K_ZS(ils,isp,ng)
                Clz(ifid)=mmol2gww*Blz*ZLpref(ils,isp,ng)/              &
     &                                    K_ZL(ils,isp,ng)
                Cpz(ifid)=mmol2gww*Bpz*ZPpref(ils,isp,ng)/              &
     &                                    K_ZP(ils,isp,ng)
                Csum=Csz(ifid)+Clz(ifid)+Cpz(ifid)
                pvalue=(Csz(ifid)+Clz(ifid)+Cpz(ifid))/(1.0_r8+Csum)
! Assimilation efficiency
                Aeff=MIN(a_AE(ils,isp,ng)*Fweight**b_AE(ils,isp,ng),    &
     &                                            AEmax(ils,isp,ng))
                fac=dtdays*Fbiomass*gww2mmol*CmaxF(ifid)*               &
     &                                Aeff/(1.0_r8+Csum)
                CszTot(k)=CszTot(k)+fac*Csz(ifid)
                ClzTot(k)=ClzTot(k)+fac*Clz(ifid)
                CpzTot(k)=CpzTot(k)+fac*Cpz(ifid)
              END IF
              thisfish => thisfish % next
            END DO
          END IF
! Check and adjust for overgrazing
          DO k=1,N(ng)
            Vcell=Hz(i,j,k)/(pm(i,j)*pn(i,j))
            IF (Fcount(k).gt.0) THEN
              BszAvg(k)=BszAvg(k)/REAL(Fcount(k),r8)
              BlzAvg(k)=BlzAvg(k)/REAL(Fcount(k),r8)
              BpzAvg(k)=BpzAvg(k)/REAL(Fcount(k),r8)
            END IF
            IF (CszTot(k).gt.(BszAvg(k)*Vcell))                         &
     &        Csz_fac(k)=BszAvg(k)*Vcell/CszTot(k)
            IF (ClzTot(k).gt.(BlzAvg(k)*Vcell))                         &
     &        Clz_fac(k)=BlzAvg(k)*Vcell/ClzTot(k)
            IF (CpzTot(k).gt.(BpzAvg(k)*Vcell))                         &
     &        Cpz_fac(k)=BpzAvg(k)*Vcell/CpzTot(k)
          END DO
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              isp = idfish(species(ifid))
              IF (fishthread(ifid).and.bounded(ifid).and.               &
     &            alive(ifid).and.                                      &
     &            (lifestage(ifid).ge.if_larva)) THEN
                Fweight=MAX(0.0_r8,bioenergy(ifwwt,ifid))
                Fworth=bioenergy(ifworth,ifid)
                Flength=bioenergy(iflngth,ifid)
                Fbiomass=Fweight*Fworth
                k=MIN(MAX(NINT(track(izgrd,itime,ifid)),1),N(ng))
                ils=lifestage(ifid)
                Csz(ifid)=Csz_fac(k)*Csz(ifid) 
                Clz(ifid)=Clz_fac(k)*Clz(ifid) 
                Cpz(ifid)=Cpz_fac(k)*Cpz(ifid) 
                Csum=Csz(ifid)+Clz(ifid)+Cpz(ifid)
                pvalue=(Csz(ifid)+Clz(ifid)+Cpz(ifid))/(1.0_r8+Csum)
                pvalue=MIN(pvalue,pvalmax(ils,isp,ng))
                Aeff=MIN(a_AE(ils,isp,ng)*Fweight**b_AE(ils,isp,ng),    &
     &                                            AEmax(ils,isp,ng))
                Csmp=CmaxF(ifid)*pvalue*Aeff
! Respiration
! Swimming velocity based on fish weight (in cm/s)
                IF (activity(ils,isp,ng).gt.0.0_r8) THEN
                  Actv=activity(ils,isp,ng)
                ELSE
                  Actv=EXP(d_R(ils,isp,ng)*Fswim(ils,isp,ng)*           &
     &                                        0.1_r8*Flength)
                END IF
                fofT4R=EXP(cr(ils,isp,ng)*(Ftemp-tr(ils,isp,ng))) 
                Resp=a_R(ils,isp,ng)*Fweight**(-b_R(ils,isp,ng))*       &
     &                                       foft4R*Actv*5.285_r8
! Fixed prey and fish density
!                preyJ=Cal_Z(isp,ng)
!                fishJ=Cal_F(isp,ng)
! Variable prey and fish density
                preyJ=Cal_Z(isp,ng)
                IF (isp.eq.if_anchovy) tday=dy_time-270.0_r8
                IF (isp.eq.if_sardine) tday=dy_time-300.0_r8
                fishJ=Cal_F(isp,ng)+1000.0_r8*                          &
     &                  COS(tday*2.0_r8*3.14159_r8/days_year)
                Cal_ZF=preyJ/fishJ
! Update Fish Weight (g wet wt)
                bioenergy(ifwwt,ifid)=MAX(0.0_r8,                       &
     &                 Fweight*(1.0_r8+dtdays*(Csmp-Resp)*Cal_ZF))
!
! Additional output variables
! Fish p-value
                bioenergy(ifpval,ifid)=pvalue
! Fish consumption
                fac=CmaxF(ifid)*Aeff/(1.0_r8+Csum)
                bioenergy(ifcsmp,ifid)=Csmp
                bioenergy(ifcsmPS,ifid)=0.0_r8
                bioenergy(ifcsmPL,ifid)=0.0_r8
                bioenergy(ifcsmZS,ifid)=fac*Csz(ifid)
                bioenergy(ifcsmZL,ifid)=fac*Clz(ifid)
                bioenergy(ifcsmZP,ifid)=fac*Cpz(ifid)
! Fish respiration
                bioenergy(ifresp,ifid)=Resp
! Compute new length
! JF NOTE: For now, J->A is at 1st birthday (1 Jan)
!                IF ((lifestage(ifid).eq.if_juvenile).and.               &
!     &              (bioenergy(ifwwt,ifid).gt.WeightJA(isp,ng)))        &
!     &              lifestage(ifid)=if_subadult
                ils=lifestage(ifid)
                LfromW=(1.0_r8/aw2l(ils,isp,ng)*                        &
     &              bioenergy(ifwwt,ifid))**(1.0_r8/bw2l(ils,isp,ng))
                deltaL1=MAX(LfromW-Flength,0.0_r8)
                deltaL2=MAX((1.0_r8-EXP(-dSLk(ils,isp,ng)*dtdays))*     &
     &                      (dSLinf(ils,isp,ng)-Flength),0.0_r8)
                bioenergy(iflngth,ifid)= Flength+MIN(deltaL1,deltaL2)
!
# ifdef FISH_FEEDBACK
! NEMURO sink terms for ZOO consumption and for DON egestion/excretion
! (will need to be converted from mmolN to mmolN/m3 in NEMURO)
! Ramp up feedback over fisrt ten years
!                t_ramp=(tdays(ng)-dstart)/(10.0_r8*days_year)
!                t_ramp=MIN(1.0_r8,t_ramp)
                t_ramp=1.0_r8
                fac=t_ramp*dtdays*Fbiomass*gww2mmol
                feedback(iSphy,ifid)=fac*bioenergy(ifcsmPS,ifid)
                feedback(iLphy,ifid)=fac*bioenergy(ifcsmPL,ifid)
                feedback(iSzoo,ifid)=fac*bioenergy(ifcsmZS,ifid)
                feedback(iLzoo,ifid)=fac*bioenergy(ifcsmZL,ifid)
                feedback(iPzoo,ifid)=fac*bioenergy(ifcsmZP,ifid)
                feedback(iPON_,ifid)=0.0_r8
# endif
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
