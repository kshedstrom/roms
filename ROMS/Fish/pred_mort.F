#include "cppdefs.h"
      MODULE pred_mort_mod
#if defined NONLINEAR && defined NEMURO_SAN && defined PREDATOR
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
!     bioenergy  Bioenergetic pred fields.                             !
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: pred_mort

      CONTAINS
!
!***********************************************************************
      SUBROUTINE pred_mort (ng, LBi, UBi, LBj, UBj, LBk, UBk, itime,    &
     &                            pm, pn,                               &
# ifdef SOLVE3D
     &                            Hz,                                   &
# endif
     &                            MyThread, bounded, track,             &
     &                            bioenergy, alive, eatfish, species,   &
     &                            pred_count, pred_list, prednodes,     &
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
      integer, intent(in) :: pred_count(LBi:UBi,LBj:UBj)
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)
      integer, intent(inout) :: species(Npred(ng))

      type(prednode), intent(in) :: pred_list(LBi:UBi,LBj:UBj)
      type(prednode), target, intent(in) :: prednodes(Npred(ng))
      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Npred(ng))

      logical, intent(in) :: MyThread(Npred(ng))
      logical, intent(inout) :: bounded(Npred(ng))
      logical, intent(inout) :: alive(Npred(ng))
      logical, intent(in) :: eatfish(Npred(ng))

      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
# endif
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Npred(ng))
      real(r8), intent(inout) :: bioenergy(NPredV(ng),Npred(ng))
!
!  Local variable declarations.
!
      integer :: i1, i2, j1, j2, i, j, ii, jj
      integer :: ipred, ipid, ifish, ifid, iFwt, ipsp, ifsp, nsp

      type(prednode), pointer :: thispred
      type(fishnode), pointer :: thisfish

      real(r8) :: Pweight, Pworth, Px, Py, Pz
      real(r8) :: Fweight, Fworth, Fx, Fy, Fz
      real(r8) :: CellVol
      real(r8) :: FTotFish
      real(r8) :: dtdays, Csmp, Egst, Excr, Resp, Actv, SDA
      real(r8) :: mmol2gww, gww2mmol, Cmax, Uswim, Fpred, dist
      real(r8) :: fac

      real(r8) :: FishVul(Nspecies(ng))
      real(r8) :: PopVul(Nspecies(ng))
      real(r8) :: BioVul(Nspecies(ng))
      real(r8) :: FFish(Nspecies(ng))
      real(r8) :: CFish(Nspecies(ng))
      real(r8) :: Pmort(Nspecies(ng))
!
!
!-----------------------------------------------------------------------
!  Compute predation mortality on fish
!-----------------------------------------------------------------------
!
      dtdays=dt(ng)*sec2day
!
      mmol2gww=1.0_r8
      gww2mmol=1.0_r8/mmol2gww
!
! Fish consumption by species from each predator
! Note: fish_by_cell and pred_by_cell have contrib from whole water column
!
      DO i=LBi+1,UBi-1
        DO j=LBj+1,UBj-1
          IF (pred_count(i,j).gt.0) THEN
            thispred => pred_list(i,j) % next
            DO ipred=1,pred_count(i,j)
              ipid = thispred % pred
              ipsp = idpred(species(ipid))
              IF (MyThread(ipid).and.bounded(ipid).and.                 &
     &                               alive(ipid)) THEN
                Pweight=bioenergy(ipwwt,ipid)
                Pworth=bioenergy(ipworth,ipid)
                Px=track(ixgrd,itime,ipid)
                Py=track(iygrd,itime,ipid)
! Swimming velocity based on fish weight (in cm/s)
                Uswim=a_Swim(ipsp,ng)*Pweight**b_Swim(ipsp,ng)

! Max. Consumption (NEEDS TEMPERATURE DEPENDENCE)
                Cmax=a_Cmax(ipsp,ng)*Pweight**(-b_Cmax(ipsp,ng))
                IF (eatfish(ipid)) THEN
                  DO nsp=1,Nspecies(ng)
                    FishVul(nsp)=0.0_r8
                    PopVul(nsp)=0.0_r8
                    BioVul(nsp)=0.0_r8
                    Pmort(nsp)=0.0_r8
                  END DO
                  DO ii=i-1,i+1
                    DO jj=j-1,j+1
                      thisfish => fish_list(ii,jj) % next
                      DO ifish=1,fish_count(ii,jj)
                        ifid = thisfish % fish
                        ifsp = idfish(species(ifid))
                        IF (MyThread(ifid).and.bounded(ifid).and.       &
     &                                         alive(ifid)) THEN
                          Fweight=bioenergy(ifwwt,ifid)
                          Fworth=bioenergy(ifworth,ifid)
                          Fx=track(ixgrd,itime,ifid)
                          Fy=track(iygrd,itime,ifid)
                          dist=(((Px-Fx)/pm(i,j))**2+                   &
     &                          ((Py-Fy)/pn(i,j))**2)**0.5_r8
                          IF (dist.lt.0.01_r8*Uswim*dt(ng)) THEN
                            FishVul(ifsp)=FishVul(ifsp)+Fweight
                            PopVul(ifsp)=PopVul(ifsp)+Fworth
                            BioVul(ifsp)=BioVul(ifsp)+Fweight*Fworth
                          END IF 
                        END IF 
                        thisfish => thisfish % next
                      END DO
                    END DO
                  END DO
! PK_fish in units of biomass (g fish) instead of concentration (g fish/m3)
                  FTotFish=0.0_r8
                  DO nsp=1,Nspecies(ng)
                    FFish(nsp)=BioVul(nsp)*                             &
     &                         Fpref(ipsp,ng)/K_Fish(ipsp,ng)
                    FTotFish=FTotFish+FFish(nsp)
                  END DO
                  DO nsp=1,Nspecies(ng)
                    CFish(nsp)=Cmax*FFish(nsp)/(1.0_r8+FTotFish)
! Must kill individual fish eaten to get correct consumption rates
! Approach: Distribute mortality evenly over vulnerable individuals in cell
                    IF (CFish(nsp).gt.0.0_r8) THEN
                      fac=CFish(nsp)*Pweight*Pworth/BioVul(nsp)
                    ELSE
                      fac=0.0_r8
                    END IF
                    IF (fac.ge.1.0_r8) THEN
                      fac=0.99_r8
                      CFish(nsp)=fac*BioVul(nsp)/(Pweight*Pworth)
                      FFish(nsp)=CFish(nsp)*                        &
     &                           (1.0_r8+FTotFish)/Cmax
                    END IF
                    Pmort(nsp)=LOG(1.0_r8-fac)
                  END DO
! Decrease worth of prey individuals accordingly
                  DO ii=i-1,i+1
                    DO jj=j-1,j+1
                      thisfish => fish_list(ii,jj) % next
                      DO ifish=1,fish_count(ii,jj)
                        ifid = thisfish % fish
                        ifsp = idfish(species(ifid))
                        IF (MyThread(ifid).and.bounded(ifid).and.       &
     &                                         alive(ifid)) THEN
                          Fweight=bioenergy(ifwwt,ifid)
                          Fworth=bioenergy(ifworth,ifid)
                          Fx=track(ixgrd,itime,ifid)
                          Fy=track(iygrd,itime,ifid)
                          dist=(((Px-Fx)/pm(i,j))**2+                   &
     &                          ((Py-Fy)/pn(i,j))**2)**0.5_r8
                          IF (dist.lt.0.01_r8*Uswim*dt(ng)) THEN
                             bioenergy(ifworth,ifid)=Fworth-            &
     &                          (1.0_r8-exp(Pmort(ifsp)*dtdays))*Fworth
! Predation mortality by species for output
!                            fmortP(isp,ofid)=fmortP(isp,ofid)+Pmort(osp)
                          END IF
                        END IF
                        thisfish => thisfish % next
                      END DO
                    END DO
                  END DO
! Recompute FTotFish based on updated FFish
                  FTotFish=0.0_r8
                  DO nsp=1,Nspecies(ng)
                    FTotFish=FTotFish+FFish(nsp)
                  END DO
                ELSE
                  DO nsp=1,Nspecies(ng)
                    FFish(nsp)=0.0_r8
                    CFish(nsp)=0.0_r8
                  END DO
                  FTotFish=0.0_r8
                END IF
!
! Pred P-value
                bioenergy(ippval,ipid)=Csmp/Cmax
! Pred consumption
                bioenergy(ipcsmF1,ipid)=CFish(1)
                bioenergy(ipcsmF2,ipid)=CFish(2)
                bioenergy(ipcsmF3,ipid)=CFish(3)
                bioenergy(ipcsmF4,ipid)=CFish(4)
                bioenergy(ipcsmF5,ipid)=CFish(5)

              END IF
              thispred => thispred % next
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE pred_mort
#endif
      END MODULE pred_mort_mod
