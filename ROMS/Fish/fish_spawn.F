#include "cppdefs.h"
      MODULE fish_spawn_mod
#ifdef NEMURO_SAN
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes numbers and locations of eggs spawned today.  !
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: fish_spawn

      CONTAINS
!
!***********************************************************************
      SUBROUTINE fish_spawn (ng, tile)
!***********************************************************************

      USE mod_param
      USE mod_types
      USE mod_fish
      USE mod_grid
      USE mod_stepping
      USE mod_biology
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
      CALL fish_spawn_tile (ng, tile, nnew(ng), LBi, UBi, LBj, UBj,     &
     &                        nfp1(ng),                                 &
     &                        GRID(ng) % spawn_dist,                    &
     &                        FISHES(ng) % bounded,                     &
     &                        FISHES(ng) % track,                       &
     &                        FISHES(ng) % bioenergy,                   &
     &                        FISHES(ng) % alive,                       &
     &                        FISHES(ng) % species,                     &
     &                        FISHES(ng) % lifestage,                   &
     &                        FISHES(ng) % num_free,                    &
     &                        FISHES(ng) % num_super)
!
      END SUBROUTINE fish_spawn
!
!***********************************************************************
      SUBROUTINE fish_spawn_tile (ng, tile, nnew, LBi, UBi, LBj, UBj,   &
     &                              nfp1, spawn_dist,                   &
     &                              bounded, track, bioenergy, alive,   &
     &                              species, lifestage,                 &
     &                              num_free, num_super)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_types
      USE mod_biology
      USE mod_parallel
      USE interp_fish_mod
      USE mod_ncparam
      USE mod_grid
      USE mod_fish
# ifdef DISTRIBUTE
      USE distribute_mod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nnew
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nfp1
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: spawn_dist(LBi:,LBj:)
      logical, intent(in) :: bounded(:)
      real(r8), intent(in) :: track(:,0:,:)
      real(r8), intent(inout) :: bioenergy(:,:)
      logical, intent(in) :: alive(:)
      integer, intent(in) :: species(:)
      integer, intent(inout) :: lifestage(:)
      integer, intent(in) :: num_free(:)
      integer, intent(out) :: num_super(:)
# else
      real(r8), intent(in) :: spawn_dist(LBi:UBi,LBj:UBj)
      logical, intent(in) :: bounded(Nfish(ng))
      real(r8), intent(in) :: track(NFV(ng),0:NFT,Nfish(ng))
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfish(ng))
      logical, intent(in) :: alive(Nfish(ng))
      integer, intent(in) :: species(Nfish(ng))
      integer, intent(inout) :: lifestage(Nfish(ng))
      integer, intent(in) :: num_free(Nspecies(ng))
      integer, intent(out) :: num_super(Nspecies(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, l, isp, ils, ifish, ifid
      integer  :: mo1, mo2

      real(r8) :: Ftemp, Fweight, Flength, Fworth, delCsmp, Csmp
      real(r8) :: WfromL, LfromW, eggs, batchJ0, extraJ, dJ
      real(r8) :: rest, ripe, rcvr, gonadJ, develop, atresia
      real(r8) :: duration, matchk, eggcostJ, eggcostWT, preyJ, fishJ
      real(r8) :: sp_time, dt1day, tday, dtspawn
      real(r8) :: dmonth, fac

      real(r8) :: Tclm(2,12)
      real(r8) :: ZSclm(2,12)
      real(r8) :: ZLclm(2,12)
      real(r8) :: ZPclm(2,12)
      real(r8) :: Ftemp_clm(2)

# ifdef DISTRIBUTE
      real(r8) :: Xstr, Xend, Ystr, Yend
      real(r8), parameter :: Fspv = 0.0_r8
      integer, parameter :: iFspv = 0
      real(r8), dimension(Nfish(ng)*NFishV(ng)) :: FwrkF
      integer,  dimension(Nfish(ng)) :: FwrkI
      integer  :: NptsF, NptsL
#endif
      logical  :: MyFishThread(Nfish(ng))
      logical, parameter :: Gmask = .FALSE.
# ifdef MASKING
      logical, parameter :: Lmask = .TRUE.
# else
      logical, parameter :: Lmask = .FALSE.
# endif

# include "set_bounds.h"
# ifdef DISTRIBUTE
!-----------------------------------------------------------------------
! In distributed-memory configuration, determine which node bounds the
! current location of the fish. Assign non-bounded fish to the
! master node.
!-----------------------------------------------------------------------
!
! The strategy here is to build a switch that processes only the fish
! contained within the node bounds. The trajectory data for the new 
! time-level (nfp1) is initialized to Fspv. These values are used during
! recombining step at the end of the routine.  Since a SUM reduction is
! carried-out, setting Fspv to zero means the fish only contribute in
! their own tile.
!
      NptsF=NFishV(ng)*Nfish(ng)
      NptsL=Nfish(ng)

      Xstr=REAL(BOUNDS(ng)%Istr(MyRank),r8)-0.5_r8
      Xend=REAL(BOUNDS(ng)%Iend(MyRank),r8)+0.5_r8
      Ystr=REAL(BOUNDS(ng)%Jstr(MyRank),r8)-0.5_r8
      Yend=REAL(BOUNDS(ng)%Jend(MyRank),r8)+0.5_r8
      DO l=1,Nfish(ng)
        MyFishThread(l)=.FALSE.
        IF ((Xstr.le.track(ixgrd,nfp1,l)).and.                          &
     &      (track(ixgrd,nfp1,l).lt.Xend).and.                          &
     &      (Ystr.le.track(iygrd,nfp1,l)).and.                          &
     &      (track(iygrd,nfp1,l).lt.Yend)) THEN
          MyFishThread(l)=.TRUE.
        ELSE IF (Master.and.(.not.bounded(l))) THEN
          MyFishThread(l)=.TRUE.
        ELSE
          DO i=1,NFishV(ng)
            bioenergy(i,l)=Fspv
          END DO
          lifestage(l)=iFspv
        END IF
      END DO
# else
      DO l=1,Nfish(ng)
        MyFishThread(l)=.TRUE.
      END DO
# endif
!
! sp_time: time in days, modulo days_year
! Note: only works if time origin is at the start of a year
      sp_time=REAL(INT(time(ng)/86400.0_r8/days_year))
      sp_time=time(ng)/86400.0_r8-days_year*sp_time
! Since the routine is call once a day, delta t must be 1 day
      dt1day=1.0_r8
!
      DO ifid=1,Nfish(ng)
        IF (bounded(ifid).and.MyFishThread(ifid).and.alive(ifid)) THEN
          isp = idfish(species(ifid))
          Flength=bioenergy(iflngth,ifid)
! Check for maturity
          IF (lifestage(ifid).eq.if_subadult) THEN
            matchk=1.0_r8/(1.0_r8+EXP(-(amature(isp,ng)+                &
     &                        bmature(isp,ng)*Flength)))
            IF (matchk.ge.0.5_r8) lifestage(ifid)=if_adult
          END IF
          ils = lifestage(ifid)
          IF (lifestage(ifid).eq.if_adult) THEN
! Variables needed for spawning
            Ftemp=track(itemp+NFV(ng)-NT(ng),nnew,ifid)
            Fweight=bioenergy(ifwwt,ifid)
            Fworth=bioenergy(ifworth,ifid)
            WfromL=al2w(ils,isp,ng)*Flength**bl2w(ils,isp,ng)
            ripe=bioenergy(ifripe,ifid)
            rest=bioenergy(ifrest,ifid)
            develop=bioenergy(ifdevl,ifid)
            atresia=bioenergy(ifatre,ifid)
            gonadJ=bioenergy(ifgonj,ifid)
            rcvr=bioenergy(ifrcvr,ifid)
            eggs=bioenergy(ifdegg,ifid)
            eggcostJ=bioenergy(ifecst,ifid)
            dJ=bioenergy(ifedjt,ifid)
            Csmp=bioenergy(ifcsmp,ifid)
! Fish energy density (varies seasonally)
            preyJ=Cal_Z(isp,ng)
            IF (isp.eq.if_anchovy) tday=sp_time-270.0_r8
            IF (isp.eq.if_sardine) tday=sp_time-300.0_r8
            fishJ=Cal_F(isp,ng)+1000.0_r8*                              &
     &              COS(tday*2.0_r8*3.14159_r8/days_year)
            IF (ripe.ge.1.0_r8) THEN
              rest=0.0_r8
              develop=0.0_r8
              ripe=0.0_r8
              eggs=0.0_r8
              gonadJ=0.0_r8
              eggcostJ=0.0_r8
              dJ=0.0_r8
            END IF
            duration=abatch(isp,ng)*EXP(-bbatch(isp,ng)*                &
     &                          (Ftemp-T0batch(isp,ng)))
            IF ((develop.le.0.1_r8).and.(atresia.le.0.1_r8)) THEN
              rest=rest-dt1day/duration
              batchJ0=epg(isp,ng)*Fweight*eegg(isp,ng)*megg(isp,ng)
              delCsmp=Csmp*(preyJ/fishJ)*Fweight*duration
              IF (INT(breed(isp,ng)+0.1).eq.1) THEN
! Capital breeder
                extraJ=MAX(0.0_r8,(Fweight-WfromL)*fishJ)
              ELSE IF (INT(breed(isp,ng)+0.1).eq.2) THEN
! Income breeder
                extraJ=MAX(0.0_r8,delCsmp*duration*fishJ)
              ELSE IF (INT(breed(isp,ng)+0.1).eq.3) THEN
! Intermediate breeder
                extraJ=MAX(0.0_r8,(Fweight-WfromL)*fishJ)+              &
     &                 MAX(0.0_r8,delCsmp*duration*fishJ)
              END IF
              IF ((sp_time.ge.Fspstr(isp,ng)).and.                      &
     &            (sp_time.le.Fspend(isp,ng)).and.                      &
     &            (Ftemp.ge.FspTmin(isp,ng)).and.                       &
     &            (Ftemp.le.FspTmax(isp,ng))) THEN
                IF ((rest.lt.0.0_r8).and.(extraJ.gt.batchJ0)) THEN
                  IF (INT(breed(isp,ng)+0.1).eq.1) THEN
! Capital breeder
                    eggs=epg(isp,ng)*Fweight
                  ELSE IF (INT(breed(isp,ng)+0.1).eq.2) THEN
! Income breeder
                    eggs=epg(isp,ng)*delCsmp
                  ELSE IF (INT(breed(isp,ng)+0.1).eq.3) THEN
! Intermediate breeder
                    eggs=epg(isp,ng)*Fweight+epg(isp,ng)*delCsmp
                  END IF
                  eggcostJ=eggs*megg(isp,ng)*eegg(isp,ng)
                  develop=1.0_r8
                  atresia=0.0_r8
                ELSE
                  eggcostJ=0.0_r8
                  gonadJ=0.0_r8
                  develop=0.0_r8
                  eggs=0.0_r8
                  ripe=0.0_r8
                  atresia=0.0_r8
                END IF
              END IF
            END IF
            IF (develop.ge.0.9_r8) THEN
              IF (Fweight.gt.(WfromL*pctxwt(isp,ng))) THEN
                eggcostJ=eggcostJ-dJ
                dJ=eggcostJ*dt1day/duration
                ripe=ripe+dt1day/duration
                gonadJ=gonadJ+dJ
                rest=0.0_r8
              ELSE
                atresia=1.0_r8
                rcvr=ripe
                develop=0.0_r8
                ripe=0.0_r8
                eggcostJ=0.0_r8
              END IF
            END IF
            IF (atresia.ge.0.9_r8) THEN
              duration=apof(isp,ng)*EXP(-bpof(isp,ng)*                  &
     &                          (Ftemp-T0pof(isp,ng)))
              rcvr=rcvr-dt1day/duration
              IF (rcvr.gt.0.0_r8) THEN
                dJ=-gonadJ*pctgain(isp,ng)
                dJ=dJ*dt1day/duration
                ripe=0.0_r8
                rest=0.0_r8
              ELSE
                dJ=-gonadJ*pctgain(isp,ng)
                dJ=dJ*dt1day/duration
                gonadJ=0.0_r8
                ripe=0.0_r8
                rest=1.0_r8+rcvr
                atresia=0.0_r8
                rcvr=0.0_r8
              END IF
            END IF
            eggcostWT=dJ/fishJ
            IF (ripe.ge.1.0_r8) THEN
              bioenergy(ifeggs,ifid)=0.5*ANINT(eggs)*MAX(0.0_r8,Fworth)
              ! 0.5 is because we assume half of fishes are female and spawn
              bioenergy(ifteggs,ifid)=bioenergy(ifteggs,ifid)+          &
     &                                  bioenergy(ifeggs,ifid)
              bioenergy(ifbatch,ifid)=bioenergy(ifbatch,ifid)+1.0_r8
              dtspawn=duration-REAL(INT(duration),r8)
            ELSE
              bioenergy(ifeggs,ifid)=0.0_r8
              dtspawn=dt1day
            END IF
            bioenergy(ifwwt,ifid)=MAX(0.0_r8,Fweight-dtspawn*eggcostWT)
! Store update spawning variables
            bioenergy(ifripe,ifid)=ripe
            bioenergy(ifrest,ifid)=rest
            bioenergy(ifdevl,ifid)=develop
            bioenergy(ifatre,ifid)=atresia
            bioenergy(ifgonj,ifid)=gonadJ
            bioenergy(ifrcvr,ifid)=rcvr
            bioenergy(ifdegg,ifid)=eggs
            bioenergy(ifedjt,ifid)=dJ
            bioenergy(ifecst,ifid)=eggcostJ
          END IF
! If past spawning season, reset attributes
          IF (sp_time.gt.Fspend(isp,ng)) THEN
            bioenergy(ifripe,ifid)=0.0_r8
            bioenergy(ifrest,ifid)=0.0_r8
            bioenergy(ifdevl,ifid)=0.0_r8
            bioenergy(ifatre,ifid)=0.0_r8
            bioenergy(ifgonj,ifid)=0.0_r8
            bioenergy(ifrcvr,ifid)=0.0_r8
            bioenergy(ifdegg,ifid)=0.0_r8
            bioenergy(ifedjt,ifid)=0.0_r8
            bioenergy(ifecst,ifid)=0.0_r8
!            bioenergy(ifeggs,ifid)=0.0_r8
!            bioenergy(ifbatch,ifid)=0.0_r8
          END IF
        END IF
      END DO
!
# ifndef EGGS_BISECTION
      CALL interp_fish (ng, LBi, UBi, LBj, UBj, 1, 1,                   &
     &                    1, Nfish(ng), nfp1, ifspwnloc, Nfish(ng),     &
     &                    r2dvar, Gmask,                                &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    spawn_dist(:,:),                              &
     &                    MyFishThread, bounded, track, bioenergy)
# endif

#ifdef DISTRIBUTE
      FwrkF=RESHAPE(bioenergy,(/NptsF/))
      CALL mp_collect (ng, iNLM, NptsF, Fspv, FwrkF)
      bioenergy=RESHAPE(FwrkF,(/NFishV(ng),Nfish(ng)/))
      CALL mp_collect_i (ng, iNLM, Nfish(ng), iFspv, lifestage)
#endif
!
      RETURN
      END SUBROUTINE fish_spawn_tile
#endif
      END MODULE fish_spawn_mod
