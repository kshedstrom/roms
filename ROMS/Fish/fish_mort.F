#include "cppdefs.h"
      MODULE fish_mort_mod
#if defined NONLINEAR && defined NEMURO_SAN
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group         Mark Hadfield   !
!    Licensed under a MIT/X style license             John M. Klinck   !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes fish mortality.                               !
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: fish_mort

      CONTAINS
!
!***********************************************************************
      SUBROUTINE fish_mort (ng, tile, LBi, UBi, LBj, UBj, LBk, UBk,     &
     &                          itime, nnew, pm, pn,                    &
     &                          fishthread, bounded, track,             &
     &                          bioenergy, alive, species,              &
     &                          lifestage, deathby,                     &
     &                          fmortN, fmortS, fmortP, fmortPsum,      &
# ifdef PREDATOR
     &                          pred_count, pred_list, prednodes,       &
     &                          pred_bioenergy, pred_species,           &
     &                          pred_track,                             &
# endif
     &                          fish_count, fish_list, fishnodes)
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
# ifdef PREDATOR
      integer, intent(in) :: pred_count(LBi:UBi,LBj:UBj)
# endif
      integer, intent(inout) :: species(Nfish(ng))
# ifdef PREDATOR
      integer, intent(inout) :: pred_species(Npred(ng))
# endif
      integer, intent(inout) :: lifestage(Nfish(ng))
      integer, intent(inout) :: deathby(Nfish(ng))

      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Nfish(ng))
# ifdef PREDATOR
      type(prednode), intent(in) :: pred_list(LBi:UBi,LBj:UBj)
      type(prednode), target, intent(in) :: prednodes(Npred(ng))
# endif

      logical, intent(in) :: fishthread(Nfish(ng))
      logical, intent(inout) :: bounded(Nfish(ng))
      logical, intent(inout) :: alive(Nfish(ng))

      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfish(ng))
# ifdef PREDATOR
      real(r8), intent(in) :: pred_track(NFV(ng),0:NFT,Npred(ng))
# endif
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfish(ng))
# ifdef PREDATOR
      real(r8), intent(in) :: pred_bioenergy(NPredV(ng),Npred(ng))
# endif
      real(r8), intent(inout) :: fmortN(Nfish(ng))
      real(r8), intent(inout) :: fmortS(Nfish(ng))
      real(r8), intent(inout) :: fmortP(Nfish(ng))
      real(r8), intent(inout) :: fmortPsum(Nfish(ng))
!
!  Local variable declarations.
!
      integer :: i1, i2, j1, j2, i, j, ii, jj, k
      integer :: ifish, ifid, isp, ifsp, ils
      integer :: ipred, ipid, ipsp, nsp

      type(fishnode), pointer :: thisfish
# ifdef PREDATOR
      type(prednode), pointer :: thispred
# endif

      real(r8) :: Fweight, Fworth, Flength, Fx, Fy
      real(r8) :: Pweight, Pworth, Px, Py, FTotFish, CmaxP
      real(r8) :: dtdays, dtsecs, dy_time, t_ramp, tday
      real(r8) :: dist, dswim, fac, LfromW, WfromL
      real(r8) :: Ndmort, Fdmort, Tdmort, Wpct

      real(r8) :: FishVul(Nspecies(ng))
      real(r8) :: PopVul(Nspecies(ng))
      real(r8) :: BioVul(Nspecies(ng))
      real(r8) :: FFish(Nspecies(ng))
      real(r8) :: CFish(Nspecies(ng))
      real(r8) :: Pmort(Nspecies(ng))

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute mortality terms for fish
!-----------------------------------------------------------------------
!
!      dtdays=dt(ng)*sec2day
!      dtsecs=4.0_r8*dt(ng)
!      dtdays=4.0_r8*dt(ng)*sec2day
      dtsecs=3600.0_r8
      dtdays=3600.0_r8*sec2day
! dy_time: time in days, modulo days_year
! Note: only works if time origin is at the start of a year
      dy_time=REAL(INT(time(ng)/86400.0_r8/days_year))
      dy_time=time(ng)/86400.0_r8-days_year*dy_time
!
! Starvation mortality
!
      DO i=Istr,Iend
        DO j=Jstr,Jend
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
                ils=lifestage(ifid)
                WfromL=al2w(ils,isp,ng)*Flength**bl2w(ils,isp,ng)
                IF (lifestage(ifid).eq.if_larva) THEN
                  Wpct=0.5_r8
                ELSE
                  Wpct=0.4_r8
                END IF
                IF (Fweight.lt.(Wpct*WfromL)) THEN
                  alive(ifid)=.FALSE.
                  IF (deathby(ifid).eq.0) deathby(ifid)=1
                END IF
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
!
# ifdef PREDATOR
! Predation mortality
!
! Zero out predation mortality from previous time step
! Zero out summed predation mortality based on output frequency
!      DO i=LBi,UBi
!        DO j=LBj,UBj
      DO i=MAX(Istr-2,1),MIN(Iend+2,Lm(ng))
        DO j=MAX(Jstr-2,1),MIN(Jend+2,Mm(ng))
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              IF (fishthread(ifid).and.bounded(ifid).and.               &
     &            alive(ifid).and.                                      &
     &            (lifestage(ifid).ge.if_subadult)) THEN
                fmortP(ifid)=0.0_r8
                IF (MOD(iic(ng)-1,nFSH(ng)).eq.0) THEN
                  fmortPsum(ifid)=0.0_r8
                END IF
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
! Predation mortality at current time step
!      DO i=LBi+1,UBi-1
!        DO j=LBj+1,UBj-1
      DO i=MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
        DO j=MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
          IF (pred_count(i,j).gt.0) THEN
            thispred => pred_list(i,j) % next
            DO ipred=1,pred_count(i,j)
              ipid = thispred % pred
              ipsp = idpred(pred_species(ipid))
! Check if predators present (April 1 to September 30)
              IF ((dy_time.ge.Pmgstr(ipsp,ng)).and.                     &
     &            (dy_time.le.Pmgend(ipsp,ng))) THEN
                Pweight=pred_bioenergy(ipwwt,ipid)
                Pworth=pred_bioenergy(ipworth,ipid)
                Px=pred_track(ixgrd,itime,ipid)
                Py=pred_track(iygrd,itime,ipid)
! Max. Consumption (NEEDS TEMPERATURE DEPENDENCE)
                CmaxP=Pcmax(ipsp,ng)
                DO nsp=1,Nspecies(ng)
                  FishVul(nsp)=0.0_r8
                  PopVul(nsp)=0.0_r8
                  BioVul(nsp)=0.0_r8
                  Pmort(nsp)=0.0_r8
                END DO
                DO ii=i-1,i+1
                  DO jj=j-1,j+1
                    IF (fish_count(ii,jj).gt.0) THEN
                      thisfish => fish_list(ii,jj) % next
                      DO ifish=1,fish_count(ii,jj)
                        ifid = thisfish % fish
                        ifsp = idfish(species(ifid))
                        IF (fishthread(ifid).and.bounded(ifid).and.     &
     &                      alive(ifid).and.                            &
     &                      (lifestage(ifid).ge.if_subadult)) THEN
                          Fweight=bioenergy(ifwwt,ifid)
                          Fworth=bioenergy(ifworth,ifid)
                          Fx=track(ixgrd,itime,ifid)
                          Fy=track(iygrd,itime,ifid)
                          dist=0.5_r8*(((Px-Fx)/pm(i,j))**2+            &
     &                          ((Py-Fy)/pn(i,j))**2)**0.5_r8
                          dswim=0.01_r8*Pswim(ipsp,ng)*dtsecs
                          IF (dist.lt.dswim) THEN
                            FishVul(ifsp)=FishVul(ifsp)+Fweight
                            PopVul(ifsp)=PopVul(ifsp)+Fworth
                            BioVul(ifsp)=BioVul(ifsp)+Fweight*Fworth
                          END IF
                        END IF
                        thisfish => thisfish % next
                      END DO
                    END IF
                  END DO
                END DO
! K_fish in units of biomass (g fish) instead of concentration (g fish/m3)
                FTotFish=0.0_r8
                DO nsp=1,Nspecies(ng)
                  FFish(nsp)=BioVul(nsp)*Fpref(nsp,ng)/K_Fish(nsp,ng)
                  FTotFish=FTotFish+FFish(nsp)
                END DO
                DO nsp=1,Nspecies(ng)
                  CFish(nsp)=CmaxP*FFish(nsp)/(1.0_r8+FTotFish)
! Must kill individual fish eaten to get correct consumption rates
! Approach: Distribute mortality evenly over vulnerable individuals in cell
! "Pmort" is the fractional biomass eaten relative to available biomass
                  IF (CFish(nsp).gt.0.0_r8) THEN
                    fac=CFish(nsp)*Pweight*Pworth/BioVul(nsp)
                  ELSE
                    fac=0.0_r8
                  END IF
                  IF (fac.ge.1.0_r8) THEN
                    fac=0.99_r8
                    CFish(nsp)=fac*BioVul(nsp)/(Pweight*Pworth)
                    FFish(nsp)=CFish(nsp)*(1.0_r8+FTotFish)/CmaxP
                  END IF
                  Pmort(nsp)=fac
                END DO
! Decrease worth of prey individuals accordingly
                DO ii=i-1,i+1
                  DO jj=j-1,j+1
                    IF (fish_count(ii,jj).gt.0) THEN
                      thisfish => fish_list(ii,jj) % next
                      DO ifish=1,fish_count(ii,jj)
                        ifid = thisfish % fish
                        ifsp = idfish(species(ifid))
                        IF (fishthread(ifid).and.bounded(ifid).and.     &
     &                      alive(ifid).and.                            &
     &                      (lifestage(ifid).ge.if_subadult)) THEN
                          Fweight=bioenergy(ifwwt,ifid)
                          Fworth=bioenergy(ifworth,ifid)
                          Fx=track(ixgrd,itime,ifid)
                          Fy=track(iygrd,itime,ifid)
                          dist=0.5_r8*(((Px-Fx)/pm(i,j))**2+            &
     &                          ((Py-Fy)/pn(i,j))**2)**0.5_r8
                          dswim=0.01_r8*Pswim(ipsp,ng)*dtsecs
                          IF (dist.lt.dswim) THEN
                             fac=1.0_r8-Pmort(ifsp)
                             bioenergy(ifworth,ifid)=Fworth*fac
! Predation mortality by species for output
                             fmortP(ifid)=fmortP(ifid)-LOG(fac)
                            IF (bioenergy(ifworth,ifid).lt.1.0_r8) THEN
                              alive(ifid)=.FALSE.
                              IF (deathby(ifid).eq.0) deathby(ifid)=3
                            END IF
                          END IF
                          fmortPsum(ifid)=fmortPsum(ifid)+fmortP(ifid)
                        END IF
                        thisfish => thisfish % next
                      END DO
                    END IF
                  END DO
                END DO
              END IF
              thispred => thispred % next
            END DO
          END IF
        END DO
      END DO
# endif
!
! Natural and fishing mortality
!
!      DO i=LBi+2,UBi-2
!        DO j=LBj+2,UBj-2
      DO i=Istr,Iend
        DO j=Jstr,Jend
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              ifid = thisfish % fish
              isp = idfish(species(ifid))
              ils = lifestage(ifid)
              IF (fishthread(ifid).and.bounded(ifid).and.               &
     &            alive(ifid).and.                                      &
     &            (lifestage(ifid).ge.if_subadult)) THEN
                Fworth=bioenergy(ifworth,ifid)
! Natural mortality rate per day
                Ndmort=Nmort(ils,isp,ng)/days_year
! Update Fish Worth
                bioenergy(ifworth,ifid)=Fworth*EXP(-Ndmort*dtdays)
! Compute yield: yield=(Fm/(Fm+Nm))*(1-EXP(-(Fm+Nm)))*Worth
!                bioenergy(ifyield,ifid)=bioenergy(ifyield,ifid)+        &
!     &            (Fdmort/Tdmort)*(1.0_r8-EXP(-Tdmort*dtdays))*Fworth
! Natural mortality
                fmortN(ifid)=Ndmort
                IF (bioenergy(ifworth,ifid).lt.1.0_r8) THEN
                  alive(ifid)=.FALSE.
                  IF (deathby(ifid).eq.0) deathby(ifid)=2
                END IF
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE fish_mort
#endif
      END MODULE fish_mort_mod
