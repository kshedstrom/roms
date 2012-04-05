#include "cppdefs.h"
      MODULE new_fish_mod
#if defined NEMURO_SAN
!
!svn $Id$
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Check to see if any fish spawned and make new superindividuals      !
!  from them.
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: new_fish

      CONTAINS
!
!***********************************************************************
      SUBROUTINE new_fish(ng)
!***********************************************************************
!
      USE mod_param
      USE mod_fish
      USE spawn_end
      USE mod_biology
      USE mod_scalars
      USE mod_stepping
      USE mod_parallel
# ifdef EGGS_TREE_FORT
      USE mod_tree
# endif
# ifdef EGGS_BISECTION
      USE mod_egglist
# endif
      USE nrutil
# ifdef DISTRIBUTE
      USE distribute_mod
# endif
      USE ran_state, ONLY: ran_seed
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
      integer :: i, isp, isp2, j, k, l, Nfound, Navail, Npts, ifsh
      integer :: count(Nspecies(ng)), count_sum
      integer :: idmaxw
      integer :: Ir, Jr
# ifdef EGGS_BISECTION
      integer, pointer :: ifish(:)
      integer, pointer :: jfish(:)
# else
      integer, pointer :: moms(:)
# endif
      real(r8), pointer :: eggs(:)
      real(r8) :: wgt
      real(r8) :: sp_time
      real(r8) :: xnudg, ynudg
      real(r8) :: maxwth, sumwth

      real(r8), parameter :: Fspv = 0.0_r8
      integer, parameter :: iFspv = 0
      logical, parameter :: lFspv = .false.
# ifdef EGGS_TREE_FORT
      type (RedBlackTree) :: tree(Nspecies(ng))
# endif
# ifdef DISTRIBUTE
      real(r8) :: Xstr, Xend, Ystr, Yend
      real(r8), dimension(Nfish(ng)*NFV(ng)*(NFT+1)) :: Fwrk
      real(r8), dimension(Nfish(ng)*NFishV(ng)) :: FwrkF
      logical,  dimension(Nfish(ng)) :: FwrkL
      integer,  dimension(Nfish(ng)) :: FwrkI
# endif
# ifdef DISTRIBUTE
      integer  :: NptsF, NptsL
#endif
      logical  :: MyThread(Nfish(ng))
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif

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
      Npts=NFV(ng)*(NFT+1)*Nfish(ng)
      NptsF=NFishV(ng)*Nfish(ng)
      NptsL=Nfish(ng)

      Xstr=REAL(BOUNDS(ng)%Istr(MyRank),r8)-0.5_r8
      Xend=REAL(BOUNDS(ng)%Iend(MyRank),r8)+0.5_r8
      Ystr=REAL(BOUNDS(ng)%Jstr(MyRank),r8)-0.5_r8
      Yend=REAL(BOUNDS(ng)%Jend(MyRank),r8)+0.5_r8
      DO l=1,Nfish(ng)
        MyThread(l)=.FALSE.
        IF (Master) THEN
          MyThread(l)=.TRUE.
        ELSE
          DO j=0,NFT
            DO i=1,NFV(ng)
              FISHES(ng) % track(i,j,l)=Fspv
            END DO
          END DO
          DO i=1,NFishV(ng)
            FISHES(ng) % bioenergy(i,l)=Fspv
          END DO
          FISHES(ng) % egg_dur(l)=Fspv
          FISHES(ng) % egg_num(l)=Fspv
          FISHES(ng) % lifestage(l)=iFspv
          FISHES(ng) % swimtype(l)=iFspv
          FISHES(ng) % species(l)=iFspv
          FISHES(ng) % cellid(l)=iFspv
          FISHES(ng) % alive(l)=lFspv
        END IF
      END DO
# else
      DO l=1,Nfish(ng)
        MyThread(l)=.TRUE.
      END DO
# endif
# ifdef DISTRIBUTE
      IF (Master) THEN
        CALL ran1 (FISHES(ng)%rwalk)
      END IF
      CALL mp_bcastf (ng, iNLM, FISHES(ng)%rwalk)
# elif defined _OPENMP
!$OMP SINGLE
      CALL ran1 (FISHES(ng)%rwalk)
!$OMP END SINGLE
# else
!     IF (Lstr.eq.1) THEN
        CALL ran1 (FISHES(ng)%rwalk)
!     END IF
# endif

!  We need to know how many new superindividuals are available for each
!  species. Only the Master gets to do any of this.
      count = 0
      IF (Master) THEN

        sp_time=REAL(INT(time(ng)/86400.0_r8/days_year))
        sp_time=time(ng)/86400.0_r8-days_year*sp_time

        DO isp=1,Nspecies(ng)

          Navail = FISHES(ng) % num_super(isp)
!         print *, "New_fish avail:", Navail
          IF ((Navail.gt.0).and.                                        &
     &        (sp_time.ge.FISHES(ng)%sp_start(isp)).and.                &
     &        (sp_time.le.FISHES(ng)%sp_end(isp))) THEN
# ifdef EGGS_BISECTION
            allocate(ifish(Navail))
            allocate(jfish(Navail))
# else
            allocate(moms(Navail))
# endif
            allocate(eggs(Navail))

!  Load up the respective data structures with the eggs
# ifdef EGGS_TREE_FORT
            call tree_init(tree(isp))

            DO l=1,Nfish(ng)
              IF (FISHES(ng) % bioenergy(ifeggs, l) .gt. 0.0_r8 .and.   &
     &             idfish(FISHES(ng) % species(l)) .eq. isp) THEN
                isp2 = idfish(FISHES(ng) % species(l))
                CALL tree_insert(tree(isp2),                            &
     &                FISHES(ng) % bioenergy(ifeggs,l),                 &
     &                FISHES(ng) % bioenergy(ifspwnloc, l), l)
                FISHES(ng) % bioenergy(ifeggs, l) = 0.0_r8
		count(isp2) = count(isp2) + 1
              END IF
            END DO
# elif defined EGGS_TREE_CXX
            CALL c_tree_init(isp)
  
            DO l=1,Nfish(ng)
              IF (FISHES(ng) % bioenergy(ifeggs, l) .gt. 0.0_r8) THEN
                isp = idfish(FISHES(ng) % species(l))
                count(isp) = count(isp) + 1
                CALL c_tree_insert(isp,                                 &
     &               FISHES(ng) % bioenergy(ifspwnloc,l),               &
     &               FISHES(ng) % bioenergy(ifeggs,l), l)
!                print *, "New_fish mom ", l,                              &
!     &                FISHES(ng) % bioenergy(ifspwnloc, l),               &
!     &                FISHES(ng) % bioenergy(ifeggs, l)
              END IF
            END DO
# elif defined EGGS_BISECTION
            CALL egglist_init(ng)
            CALL egglist_split(ng, isp, Navail, count(isp))
# endif

!             call tree_traverse(tree(is))
            IF (count(isp) > 0) THEN
# ifdef EGGS_TREE_FORT
	      print *, 'New_fish found ', count(isp),                   &
     &             'spawners in species ', isp
              CALL tree_collect(tree(isp), Navail,                      &
     &             Nfound, eggs, moms)
# elif defined EGGS_TREE_CXX
              CALL c_tree_collect(isp, FISHES(ng) % num_super(isp),     &
     &             Nfound, eggs, moms)
# elif defined EGGS_BISECTION
              CALL egglist_collect(isp, FISHES(ng) % num_super(isp),    &
     &             Nfound, eggs, ifish, jfish)
# endif
!             do i=1,Nfound
!               print *, i, eggs(i), moms(i)
!             end do
              DO i=1,Nfound
                ifsh = FISHES(ng) % next_free(isp)
                FISHES(ng) % next_free(isp) =                           &
     &                  FISHES(ng) % next_free(isp) + 1
                FISHES(ng) % num_free(isp) =                            &
     &                  FISHES(ng) % num_free(isp) - 1
                FISHES(ng) % bounded(ifsh) = .TRUE.
!  Place it at the mother's location
                xnudg = FISHES(ng) % rwalk(ifsh)
                ynudg = FISHES(ng) % rwalk(ifsh + Nfish(ng))
! JF: NEED TO MAKE SURE I,J IS NOT A LAND CELL !!!!!!!!
                DO j=0,NFT
# ifdef EGGS_BISECTION
                  FISHES(ng)%track(ixgrd,j,ifsh)=ifish(i)+              &
     &                                           (xnudg-0.5_r8)
                  FISHES(ng)%track(iygrd,j,ifsh)=jfish(i)+              &
     &                                           (ynudg-0.5_r8)
                  FISHES(ng)%track(izgrd,j,ifsh)=FISHES(ng)%Fz0(ifsh)
# else
!                  DO k=1,NFV(ng)
!                    FISHES(ng) % track(k,j,ifsh) =                      &
!     &                 FISHES(ng) % track(k,j,moms(i))
!                  END DO
                  FISHES(ng)%track(ixgrd,j,ifsh) =                      &
     &                       FISHES(ng)%track(ixgrd,j,moms(i))
                  FISHES(ng)%track(iygrd,j,ifsh) =                      &
     &                       FISHES(ng)%track(iygrd,j,moms(i))
                  FISHES(ng)%track(izgrd,j,ifsh) =                      &
     &                       FISHES(ng)%track(izgrd,j,moms(i))
# endif
                  FISHES(ng)%track(ixrhs,j,ifsh)=0.0_r8
                  FISHES(ng)%track(iyrhs,j,ifsh)=0.0_r8
                  FISHES(ng)%track(izrhs,j,ifsh)=0.0_r8
                END DO
# ifdef EGGS_BISECTION
                FISHES(ng)%Tinfo(ixgrd,ifsh)=ifish(i)+(xnudg-0.5_r8)
                FISHES(ng)%Tinfo(iygrd,ifsh)=jfish(i)+(ynudg-0.5_r8)
                FISHES(ng)%Tinfo(izgrd,ifsh)=FISHES(ng)%Fz0(ifsh)
# else
                FISHES(ng)%Tinfo(ixgrd,ifsh) =                          &
     &                     FISHES(ng)%Tinfo(ixgrd,moms(i))
                FISHES(ng)%Tinfo(iygrd,ifsh) =                          &
     &                     FISHES(ng)%Tinfo(iygrd,moms(i))
                FISHES(ng)%Tinfo(izgrd,ifsh) =                          &
     &                     FISHES(ng)%Tinfo(izgrd,moms(i))
# endif
! age zero eggs...
                FISHES(ng) % bioenergy(ifworth,ifsh) = eggs(i)
                FISHES(ng) % bioenergy(ifwwt,ifsh) = 0.0_r8
                FISHES(ng) % bioenergy(ifage,ifsh) = 0.0_r8
                FISHES(ng) % species(ifsh) = idfish_inv(isp)
                FISHES(ng) % alive(ifsh) = .TRUE.
                FISHES(ng) % lifestage(ifsh) = if_egg
                FISHES(ng) % swimtype(ifsh) = Fbehave(isp,ng)
! number individuals entering egg stage (for life table)
                FISHES(ng)%egg_num(ifsh)=eggs(i)
                FISHES(ng)%egg_dur(ifsh)=0.0_r8
! Save initial worth for splitting at end of spawning season
                FISHES(ng) % bioenergy(ifiniwth,ifsh) =                 &
     &                       FISHES(ng) % bioenergy(ifworth,ifsh)
! Cell identification for each fish (needed for fishing fleet)
                Ir=NINT(FISHES(ng)%Tinfo(ixgrd,ifsh))
                Jr=NINT(FISHES(ng)%Tinfo(iygrd,ifsh))
                FISHES(ng)%cellid(ifsh)=Ir+(Lm(ng)+2)*Jr
!
# ifdef EGGS_BISECTION
                print *, 'New_fish spawn: ', ifsh, eggs(i),             &
     &                                 ifish(i)+xnudg, jfish(i)+ynudg
# else
                print *, 'New_fish spawn: ', ifsh, eggs(i), moms(i)
# endif
                ifsh = ifsh + 1
              END DO
            END IF
# ifdef EGGS_TREE_FORT
            CALL tree_destroy(tree(1))
# elif defined EGGS_TREE_CXX
            CALL c_tree_trim
# elif defined EGGS_BISECTION
            CALL egglist_destroy(1)
# endif
# ifdef EGGS_BISECTION
            deallocate(ifish)
            deallocate(jfish)
# else
            deallocate(moms)
# endif
            deallocate(eggs)
          END IF
! Allocate remaining individuals at the end of the spawning season
! based on worth at creation
          IF ((FISHES(ng)%num_free(isp).gt.0).and.                      &
     &        (FISHES(ng)%num_free(isp).lt.Nfishperyear(ng)).and.       &
     &        (sp_time.gt.FISHES(ng)%sp_end(isp))) THEN
            print*, 'SPLIT', FISHES(ng)%num_free(isp)
            CALL split_biggest (ng, isp, FISHES(ng)%num_free(isp))
          END IF
        END DO
      END IF

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect fish on all nodes.
!-----------------------------------------------------------------------
!
      Fwrk=RESHAPE(FISHES(ng) % track,(/Npts/))
      CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
      FISHES(ng) % track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Nfish(ng)/))

      FwrkF=RESHAPE(FISHES(ng) % bioenergy,(/NptsF/))
      CALL mp_collect (ng, iNLM, NptsF, Fspv, FwrkF)
      FISHES(ng) % bioenergy=RESHAPE(FwrkF,(/NFishV(ng),Nfish(ng)/))
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv,                       &
     &                                   FISHES(ng) % egg_dur)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv,                       &
     &                                   FISHES(ng) % egg_num)
      CALL mp_collecti (ng, iNLM, Nfish(ng), iFspv,                     &
     &                                   FISHES(ng) % lifestage)
      CALL mp_collecti (ng, iNLM, Nfish(ng), iFspv,                     &
     &                                   FISHES(ng) % species)
      CALL mp_collecti (ng, iNLM, Nfish(ng), iFspv,                     &
     &                                   FISHES(ng) % swimtype)
      CALL mp_collecti (ng, iNLM, Nfish(ng), iFspv,                     &
     &                                   FISHES(ng) % cellid)
      CALL mp_collectl (ng, iNLM, Nfish(ng), FISHES(ng) % alive)
      CALL mp_collectl (ng, iNLM, Nfish(ng), FISHES(ng) % bounded)
# endif

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE new_fish
!
#endif
      END MODULE new_fish_mod
