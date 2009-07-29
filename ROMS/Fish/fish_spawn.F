#include "cppdefs.h"
      MODULE fish_spawn_mod
#ifdef NEMURO_SAN
!
!svn $Id: interp_floats.F 294 2009-01-09 21:37:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group         Mark Hadfield   !
!    Licensed under a MIT/X style license             John M. Klinck   !
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
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
      CALL fish_spawn_tile (ng, tile, nnew(ng), LBi, UBi, LBj, UBj,     &
     &                        nf(ng),                                   &
     &                        GRID(ng) % spawn_dist,                    &
     &                        FISHES(ng) % bounded,                     &
     &                        FISHES(ng) % track,                       &
     &                        FISHES(ng) % bioenergy,                   &
     &                        FISHES(ng) % alive,                       &
     &                        FISHES(ng) % species,                     &
     &                        FISHES(ng) % sp_start,                    &
     &                        FISHES(ng) % sp_end,                      &
     &                        FISHES(ng) % num_free,                    &
     &                        FISHES(ng) % num_super)
!
      END SUBROUTINE fish_spawn
!
!***********************************************************************
      SUBROUTINE fish_spawn_tile (ng, tile, nnew, LBi, UBi, LBj, UBj,   &
     &                              nf, spawn_dist,                     &
     &                              bounded, track, bioenergy, alive,   &
     &                              species, sp_start, sp_end,          &
     &                              num_free, num_super)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_types
      USE mod_biology
      USE mod_parallel
# ifdef DISTRIBUTE
      USE distribute_mod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nnew
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nf
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: spawn_dist(LBi:,LBj:)
      logical, intent(in) :: bounded(:)
      real(r8), intent(in) :: track(:,0:,:)
      real(r8), intent(inout) :: bioenergy(:,:)
      logical, intent(in) :: alive(:)
      integer, intent(in) :: species(:)
      real(r8), intent(in) :: sp_start(:)
      real(r8), intent(in) :: sp_end(:)
      integer, intent(in) :: num_free(:)
      integer, intent(out) :: num_super(:)
# else
      real(r8), intent(in) :: spawn_dist(LBi:UBi,LBj:UBj)
      logical, intent(in) :: bounded(Nfish(ng))
      real(r8), intent(in) :: track(NFV(ng),0:NFT,Nfish(ng))
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfish(ng))
      logical, intent(in) :: alive(Nfish(ng))
      integer, intent(in) :: species(Nfish(ng))
      real(r8), intent(in) :: sp_start(Nspecies(ng))
      real(r8), intent(in) :: sp_end(Nspecies(ng))
      integer, intent(in) :: num_free(Nspecies(ng))
      integer, intent(out) :: num_super(Nspecies(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, l, isp, ifish, ifid

      real(r8) :: Ftemp, temp_opt, temp_sigma
      real(r8) :: sp_time
# ifdef DISTRIBUTE
      real(r8) :: Xstr, Xend, Ystr, Yend
      real(r8), parameter :: Fspv = 0.0_r8
      real(r8), dimension(Nfish(ng)*NFishV(ng)) :: FwrkF
      integer  :: NptsF
#endif
      logical  :: MyThread(Nfish(ng))

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

      Xstr=REAL(BOUNDS(ng)%Istr(MyRank),r8)-0.5_r8
      Xend=REAL(BOUNDS(ng)%Iend(MyRank),r8)+0.5_r8
      Ystr=REAL(BOUNDS(ng)%Jstr(MyRank),r8)-0.5_r8
      Yend=REAL(BOUNDS(ng)%Jend(MyRank),r8)+0.5_r8
      DO l=1,Nfish(ng)
        MyThread(l)=.FALSE.
        IF ((Xstr.le.track(ixgrd,nf,l)).and.                            &
     &      (track(ixgrd,nf,l).lt.Xend).and.                            &
     &      (Ystr.le.track(iygrd,nf,l)).and.                            &
     &      (track(iygrd,nf,l).lt.Yend)) THEN
          MyThread(l)=.TRUE.
        ELSE IF (Master.and.(.not.bounded(l))) THEN
          MyThread(l)=.TRUE.
        ELSE
          DO i=1,NFishV(ng)
            bioenergy(i,l)=Fspv
          END DO
        END IF
      END DO
# else
      DO l=1,Nfish(ng)
        MyThread(l)=.TRUE.
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute feeding terms on zooplankton to update growth
!-----------------------------------------------------------------------
!
      temp_opt=15.0_r8
      temp_sigma=3.0_r8
!
! sp_time: time in days, modulo 360 days
      sp_time=REAL(INT(time(ng)/86400.0_r8/360.0_r8))
      sp_time=time(ng)/86400.0_r8-360.0_r8*sp_time
!
! find out how many free superindividuals per species for today
      DO isp=1, Nspecies(ng)
        num_super(isp) = num_free(isp) / (sp_end(isp) - sp_time + 1)
      END DO
!
      DO ifid=1,Nfish(ng)
        IF (bounded(ifid).and.alive(ifid).and.MyThread(ifid)) THEN
          isp = idfish(species(ifid))
! If inside spawning window compute spawning attributes
          IF ((sp_time.ge.sp_start(isp)).and.                           &
     &              (sp_time.le.sp_end(isp))) THEN
            Ftemp=track(itemp+NFV(ng)-NT(ng),nnew,ifid)
            IF (alive(ifid)) THEN
!           IF ((alive(ifid)).and.                                      &
!    &              (bioenergy(ifage,ifid).ge.Fage0(ng)).and.           &
!    &              (INT(bioenergy(ifbatch,ifid)).lt.Nbatch(ng)).and.   &
!    &              (ABS(Ftemp-temp_opt).le.temp_sigma).and.            &
!    &              ((sp_time-bioenergy(iftspwn,ifid)).ge.7.0_r8)) THEN
              bioenergy(ifeggs,ifid)=20000.0_r8*                  &
     &                                     bioenergy(ifworth,ifid)
              bioenergy(ifspwnloc,ifid)=spawn_dist(i,j)
              bioenergy(ifbatch,ifid)=bioenergy(ifbatch,ifid)+    &
     &                                      1.0_r8
              bioenergy(iftspwn,ifid)=sp_time
            ELSE
              bioenergy(ifeggs,ifid)=0.0_r8
              bioenergy(ifspwnloc,ifid)=0.0_r8
            END IF
          END IF
! If after end of spawning window, zero out spawning attributes
          IF (sp_time.gt.sp_end(isp)) THEN
            bioenergy(ifeggs,ifid)=0.0_r8
            bioenergy(ifspwnloc,ifid)=0.0_r8
          END IF
!
        END IF
      END DO
#ifdef DISTRIBUTE
      FwrkF=RESHAPE(bioenergy,(/NptsF/))
      CALL mp_collect (ng, iNLM, NptsF, Fspv, FwrkF)
      bioenergy=RESHAPE(FwrkF,(/NFishV(ng),Nfish(ng)/))
#endif
!
      RETURN
      END SUBROUTINE fish_spawn_tile
#endif
      END MODULE fish_spawn_mod
