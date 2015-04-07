#include "cppdefs.h"
      MODULE new_year_mod
#if defined NEMURO_SAN
!
!svn $Id$
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  At the end of each year, we kill off any remaining fish of the
!  the oldest age class to make room for the next age class. We also
!  move the indices for creating new fish.
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: new_year

      CONTAINS
!
!***********************************************************************
      SUBROUTINE new_year(ng, Lstr, Lend)
!***********************************************************************
!
      USE mod_param
      USE mod_fish
      USE mod_scalars
      USE mod_stepping
      USE mod_parallel
      USE mod_biology
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_collect, mp_collect_i, mp_collect_l
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend
      integer :: i, isp, isp2, j, l, Nfound, Navail, ifsh, ifid
      integer :: count(Nspecies(ng)), count_sum
      logical, dimension(Lstr:Lend) :: MyFishThread
# ifdef DISTRIBUTE
      real(r8) :: Xstr, Xend, Ystr, Yend
      real(r8), dimension(Nfish(ng)*NFV(ng)*(NFT+1)) :: Fwrk 
      real(r8), dimension(Nfish(ng)*NFishV(ng)) :: FwrkF
      real(r8), dimension(Nfish(ng)*Nspecies(ng)) :: FwrkM
      logical,  dimension(Nfish(ng)) :: FwrkL
      integer,  dimension(Nfish(ng)) :: FwrkI

      real(r8), parameter :: Fspv = 0.0_r8
      integer, parameter :: iFspv = 0
      logical, parameter :: lFspv = .false.
      integer :: NptsF
# endif
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif

# ifdef DISTRIBUTE
!
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
      IF (Lfish(ng)) THEN
        DO l=Lstr,Lend
          MyFishThread(l)=.FALSE.
          IF ((Xstr.le.FISHES(ng)%track(ixgrd,nf(ng),l)).and.           &
     &        (FISHES(ng)%track(ixgrd,nf(ng),l).lt.Xend).and.           &
     &        (Ystr.le.FISHES(ng)%track(iygrd,nf(ng),l)).and.           &
     &        (FISHES(ng)%track(iygrd,nf(ng),l).lt.Yend)) THEN
            MyFishThread(l)=.TRUE.
          ELSE IF (Master.and.(.not.FISHES(ng)%bounded(l))) THEN
            MyFishThread(l)=.TRUE.
          ELSE
            DO i=1,NFishV(ng)
              FISHES(ng)%bioenergy(i,l)=Fspv
            END DO
            FISHES(ng)%lifestage(l)=iFspv
            FISHES(ng)%deathby(l)=iFspv
            FISHES(ng)%suba_num(l)=Fspv
            FISHES(ng)%larv_num(l)=Fspv
            FISHES(ng)%larv_dur(l)=Fspv
            FISHES(ng)%alive(l)=lFspv
            MyFishThread(l)=.FALSE.
          END IF
        END DO
      END IF
# else
      DO l=1,NFish(ng)
        MyFishThread(l)=.TRUE.
      END DO
# endif
!

!  Update the indices
      DO isp=1,Nspecies(ng)
        FISHES(ng)%age_base(isp)=FISHES(ng)%age_base(isp)+              &
     &                                        Nfishperyear(ng)
        IF (FISHES(ng)%age_base(isp).ge.FISHES(ng)%species_base(isp)+   &
     &                            Nfishperyear(ng)*Nyearclass(ng)) THEN
          FISHES(ng)%age_base(isp)=FISHES(ng)%species_base(isp)
        END IF
      END DO

!  Killing the old fish for each species.
!  The max age is the same for all species (for now).

      DO isp=1,Nspecies(ng)
        DO i=1,Nfishperyear(ng)
          ifsh = i+FISHES(ng) % age_base(isp)
          FISHES(ng) % alive(ifsh) = .FALSE.
          FISHES(ng) % bounded(ifsh) = .FALSE.
          FISHES(ng) % bioenergy(ifwwt,ifsh) = 0.0_r8
          FISHES(ng) % bioenergy(iflngth,ifsh) = 0.0_r8
          FISHES(ng) % bioenergy(ifworth,ifsh) = 0.0_r8
          FISHES(ng) % bioenergy(ifiniwth,ifsh) = 0.0_r8
          FISHES(ng) % bioenergy(ifage,ifsh) = 0.0_r8
          FISHES(ng) % bioenergy(ifbday,ifsh) = 0.0_r8
          FISHES(ng) % lifestage(ifsh) = 0
          FISHES(ng) % deathby(ifsh) = 0
          DO j=0,NFT
            FISHES(ng) % track(ixgrd,j,ifsh) = 0.
            FISHES(ng) % track(iygrd,j,ifsh) = 0.
            FISHES(ng) % track(izgrd,j,ifsh) = 0.
          END DO
        END DO
      END DO

!  Get ready for the new eggs
      DO isp=1,Nspecies(ng)
        FISHES(ng)%num_free(isp)=Nfishperyear(ng)
        FISHES(ng)%num_super(isp) = FISHES(ng)%num_free(isp) /          &
     &        (Fspend(isp,ng) - Fspstr(isp,ng) + 1)
        FISHES(ng)%next_free(isp)=1+FISHES(ng)%age_base(isp)
      END DO
      IF (Master) THEN
        print *, 'New Year for Fish'
        print *, 'Free fish ', FISHES(ng) % num_free
        print *, 'Fish ages ', FISHES(ng) % age_base
      END IF

! Update fish birthday
      DO ifid=Lstr,Lend
        IF (MyFishThread(ifid).and.FISHES(ng)%bounded(ifid).and.        &
     &      FISHES(ng) % alive(ifid)) THEN
          FISHES(ng) % bioenergy(ifbday,ifid) =                         &
     &               FISHES(ng) % bioenergy(ifbday,ifid)+1.0_r8
! reset spawning variables
          FISHES(ng) % bioenergy(ifeggs,ifid)=0.0_r8
          FISHES(ng) % bioenergy(ifteggs,ifid)=0.0_r8
          FISHES(ng) % bioenergy(ifbatch,ifid)=0.0_r8
! upgrade juvenile to subadult upon reaching 1st birthday
! and record number of juvenile entering subadult stage
          IF ((FISHES(ng)%lifestage(ifid).eq.if_juvenile).and.          &
     &        (FISHES(ng)%bioenergy(ifbday,ifid).ge.1.0_r8)) THEN
            FISHES(ng) % suba_num(ifid)=                                &
     &                   FISHES(ng) % bioenergy(ifworth,ifid)
            FISHES(ng) % lifestage(ifid)=if_subadult
          END IF
! kill off larvae that did not become juveniles by end of 1st year
          IF ((FISHES(ng)%lifestage(ifid).eq.if_larva).and.             &
     &        (FISHES(ng)%bioenergy(ifbday,ifid).ge.1.0_r8)) THEN
!            FISHES(ng) % lifestage(ifid)=0
            FISHES(ng) % alive(ifid)=.FALSE.
            FISHES(ng) % bounded(ifid)=.FALSE.
            FISHES(ng) % deathby(ifid)=5
!            FISHES(ng) % bioenergy(ifwwt,ifid)=0.0_r8
!            FISHES(ng) % bioenergy(ifworth,ifid)=0.0_r8
!            FISHES(ng) % bioenergy(ifage,ifid)=0.0_r8
!            FISHES(ng) % bioenergy(ifbday,ifid)=0.0_r8
            FISHES(ng) % larv_num(ifid)=0.0_r8
            FISHES(ng) % larv_dur(ifid)=0.0_r8
          END IF
        END IF
      END DO

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect fish on all nodes.
!-----------------------------------------------------------------------
!
      FwrkF=RESHAPE(FISHES(ng)%bioenergy,(/NptsF/))
      CALL mp_collect (ng, iNLM, NptsF, Fspv, FwrkF)
      FISHES(ng)%bioenergy=RESHAPE(FwrkF,(/NFishV(ng),Nfish(ng)/))
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, FISHES(ng)%suba_num)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, FISHES(ng)%larv_dur)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, FISHES(ng)%larv_num)
      CALL mp_collect_i(ng, iNLM, Nfish(ng), iFspv,                     &
     &                 FISHES(ng)%lifestage)
      CALL mp_collect_i(ng, iNLM, Nfish(ng), iFspv, FISHES(ng)%deathby)
      CALL mp_collect_l (ng, iNLM, Nfish(ng), FISHES(ng)%alive)
!  Collect the bounded status switch.
      Fwrk=Fspv
      DO l=1,Nfish(ng)
        IF (FISHES(ng) % bounded(l)) THEN
          Fwrk(l)=1.0_r8
        END IF
      END DO
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, Fwrk)
      DO l=1,Nfish(ng)
        IF (Fwrk(l).ne.Fspv) THEN
          FISHES(ng) % bounded(l)=.TRUE.
        ELSE
          FISHES(ng) % bounded(l)=.FALSE.
        END IF
      END DO
!
# endif

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE new_year
!
#endif
      END MODULE new_year_mod
