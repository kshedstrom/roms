#include "cppdefs.h"
      MODULE new_fish_mod
#if defined NEMURO_SAN
!
!svn $Id: step_floats.F 997 2009-06-04 17:42:20Z kate $
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
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
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
      integer :: i, isp, isp2, j, k, l, Nfound, Navail, Npts, ifsh
      integer :: count(Nspecies(ng)), count_sum
      integer, pointer :: moms(:)
      real(r8), pointer :: eggs(:)
      real(r8) :: wgt

      real(r8), parameter :: Fspv = 0.0_r8
      integer, parameter :: iFspv = 0
      logical, parameter :: lFspv = .false.
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
          FISHES(ng) % lifestage(l)=iFspv
          FISHES(ng) % swimtype(l)=iFspv
          FISHES(ng) % species(l)=iFspv
          FISHES(ng) % alive(l)=lFspv
        END IF
      END DO
# else
      DO l=1,Nfish(ng)
        MyThread(l)=.TRUE.
      END DO
# endif

!  We need to know how many new superindividuals are available for each
!  species. Only the Master gets to do any of this.
      count = 0
      IF (Master) THEN

        Navail = sum(FISHES(ng) % num_super)
        IF (Navail > 0) THEN
          allocate(moms(Navail))
          allocate(eggs(Navail))

# ifdef EGGS_TREE_FORT
!  This option is currently broken...
          DO isp2=1,Nspecies(ng)
            call tree_init

            DO l=1,Nfish(ng)
              isp = idfish(FISHES(ng) % species(l))
              IF (isp2 .eq. isp .and.                                   &
     &              FISHES(ng) % bioenergy(ifeggs, l) .ne. 0) THEN
                CALL tree_insert(FISHES(ng) % bioenergy(ifspwnloc,l),   &
     &                FISHES(ng) % bioenergy(ifeggs, l), l)
                FISHES(ng) % bioenergy(ifeggs, l) = 0.0_r8
              END IF
            END DO

!           call tree_traverse
            CALL tree_collect(FISHES % num_super(isp),                  &
     &                   Nfound, eggs, moms)
            print *, 'after collect: ', Nfound
            do i=1,Nfound
              print *, i, eggs(i), moms(i)
            end do
            CALL tree_destroy
          END DO
# elif defined EGGS_TREE_CXX
          DO isp=1,Nspecies(ng)
            CALL c_tree_init(isp)
          END DO

          DO l=1,Nfish(ng)
            IF (FISHES(ng) % bioenergy(ifeggs, l) .ne. 0) THEN
              isp = idfish(FISHES(ng) % species(l))
              count(isp) = count(isp) + 1
              CALL c_tree_insert(isp,                                   &
     &             FISHES(ng) % bioenergy(ifspwnloc,l),                 &
     &             FISHES(ng) % bioenergy(ifeggs,l), l)
!              print *, "New_fish mom ", l,                              &
!     &              FISHES(ng) % bioenergy(ifspwnloc, l),               &
!     &              FISHES(ng) % bioenergy(ifeggs, l)
            END IF
          END DO

          DO isp=1,Nspecies(ng)
!           CALL c_tree_traverse(isp)
            IF (count(isp) > 0) THEN
              CALL c_tree_collect(isp, FISHES(NG) % num_super(isp),     &
     &              Nfound, eggs, moms)
              print *, 'New_fish count: ', Nfound
              DO i=1,Nfound
                ifsh = FISHES(ng) % next_free(isp)
		FISHES(ng) % next_free(isp) =                           &
     &                  FISHES(ng) % next_free(isp) + 1
                FISHES(ng) % bounded(ifsh) = .TRUE.
!  Place it at the mother's location
                DO j=0,NFT
		  DO k=1,NFV(ng)
                    FISHES(ng) % track(k,j,ifsh) =                      &
     &                 FISHES(ng) % track(k,j,moms(i))
                  END DO
                END DO
! age zero eggs...
                FISHES(ng) % bioenergy(ifworth,ifsh) = eggs(i)
                FISHES(ng) % bioenergy(ifwwt,ifsh) = 0.0_r8
                FISHES(ng) % bioenergy(ifage,ifsh) = 0.0_r8
                FISHES(ng) % species(ifsh) = idfish_inv(isp)
                FISHES(ng) % alive(ifsh) = .TRUE.
                FISHES(ng) % lifestage(ifsh) = if_egg
                FISHES(ng) % eatfish(ifsh) = .FALSE.
                FISHES(ng) % swimtype(ifsh) = 1   ! 1=Fitness, 2=Kinesis
                print *, 'New_fish: ', ifsh, eggs(i), moms(i)
                FISHES(ng) % num_free(isp) =                            &
     &                        FISHES(ng) % num_free(isp) - 1
                ifsh = ifsh + 1
              END DO
	    END IF
          END DO
          CALL c_tree_trim
# elif defined EGGS_BISECTION
          DO isp=1,Nspecies(ng)
            CALL egglist_init
            CALL egglist_split(FISHES % num_super(isp),                 &
     &                   Nfound, eggs, moms)
            CALL egglist_destroy
          END DO
# endif
          deallocate(moms)
          deallocate(eggs)
        END IF
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
      CALL mp_collecti (ng, iNLM, Nfish(ng), iFspv,                     &
     &                                   FISHES(ng) % lifestage)
      CALL mp_collecti (ng, iNLM, Nfish(ng), iFspv,                     &
     &                                   FISHES(ng) % species)
      CALL mp_collecti (ng, iNLM, Nfish(ng), iFspv,                     &
     &                                   FISHES(ng) % swimtype)
      CALL mp_collectl (ng, iNLM, Nfish(ng), FISHES(ng) % alive)
      CALL mp_collectl (ng, iNLM, Nfish(ng), FISHES(ng) % bounded)
      CALL mp_collectl (ng, iNLM, Nfish(ng), FISHES(ng) % eatfish)
# endif

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE new_fish
!
#endif
      END MODULE new_fish_mod
