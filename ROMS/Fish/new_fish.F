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
      USE nrutil
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_collect
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
      integer :: i, isp, j, l, Nfound, Navail, Npts, ifsh, count
      integer, pointer :: moms(:)
      real(r8), pointer :: eggs(:)
      real(r8) :: wgt
# ifdef DISTRIBUTE
      real(r8) :: Xstr, Xend, Ystr, Yend
      real(r8), dimension(Nfloats(ng)*NFV(ng)*(NFT+1)) :: Fwrk
      real(r8), parameter :: Fspv = 0.0_r8
# endif
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif

!  We need to know how many new superindividuals are available for each
!  species.
      count = 0
      DO l=1,Nfish(ng)
        IF (FISHES(ng) % bioenergy(ifeggs, l) .ne. 0) THEN
          count = count + 1
        END IF
      END DO
      IF (count .eq. 0) THEN
        RETURN
      END IF
      print *, 'In New_fish - found some moms with eggs:', count

      Navail = sum(FISHES(ng) % num_super)
      IF (Navail > 0) THEN
        IF (Master) THEN
	  print *, "In new_fish guts...", FISHES(ng) % num_super
          allocate(moms(Navail))
          allocate(eggs(Navail))

# ifdef EGGS_TREE_FORT
!  This option is currently broken...
          call tree_init

          DO l=1,Nfish(ng)
            IF (FISHES(ng) % bioenergy(ifeggs, l) .ne. 0) THEN
              CALL tree_insert(FISHES(ng) % bioenergy(ifspwnloc,l),     &
     &                FISHES(ng) % bioenergy(ifeggs, l), l)
            END IF
          END DO

          call tree_traverse
          CALL tree_collect(num_super, Nfound, eggs, moms)
          print *, 'after collect: ', Nfound
          do i=1,Nfound
            print *, i, eggs(i), moms(i)
          end do
          CALL tree_destroy
        END IF
# elif defined EGGS_TREE_CXX
          DO isp=1,Nspecies(ng)
            CALL c_tree_init(isp)
          END DO

          DO l=1,Nfish(ng)
            IF (FISHES(ng) % bioenergy(ifeggs, l) .ne. 0) THEN
              isp = FISHES(ng) % species(l)
              CALL c_tree_insert(isp,                                   &
     &             FISHES(ng) % bioenergy(ifspwnloc,l),                 &
     &             FISHES(ng) % bioenergy(ifeggs,l), l)
            END IF
          END DO

          DO isp=1,Nspecies(ng)
            CALL c_tree_traverse(isp)
            CALL c_tree_collect(isp, FISHES(NG) % num_super(isp),       &
     &              Nfound, eggs, moms)
            print *, 'Making ', Nfound, ' new fish of species ', isp
            ifsh = FISHES(ng) % next_free(isp)
            DO i=1,Nfound
              FISHES(ng) % bounded(ifsh) = .TRUE.
!  Place it at the mother's location
              DO j=0,NFT
                FISHES(ng) % track(ixgrd,j,ifsh) =                      &
     &                 FISHES(ng) % track(ixgrd,j,moms(i))
                FISHES(ng) % track(iygrd,j,ifsh) =                      &
     &                 FISHES(ng) % track(ixgrd,j,moms(i))
#  ifdef SOLVE3D
                FISHES(ng) % track(izgrd,j,ifsh) =                      &
     &                 FISHES(ng) % track(ixgrd,j,moms(i))
#  endif
              END DO
! age zero eggs...
              FISHES(ng) % bioenergy(ifworth,ifsh) = eggs(i)
              CALL ran1(wgt)
              FISHES(ng) % bioenergy(ifwwt,ifsh) = 2._r8*(wgt+0.5_r8)
              FISHES(ng) % bioenergy(ifage,ifsh) = 0._r8
              FISHES(ng) % species(ifsh) = isp
              FISHES(ng) % alive(ifsh) = .TRUE.
              FISHES(ng) % lifestage(ifsh) = if_egg
              print *, i, eggs(i), moms(i)
              FISHES(ng) % num_free(isp)=FISHES(ng) % num_free(isp) - 1
              ifsh = ifsh + 1
            END DO
          END DO
          CALL c_tree_destroy
# endif
          deallocate(moms)
          deallocate(eggs)
        END IF

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect fish on all nodes.
!-----------------------------------------------------------------------
!
        Npts=NFV(ng)*(NFT+1)*Nfloats(ng)
        Fwrk=RESHAPE(FISHES(ng) % track,(/Npts/))
        CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
        FISHES(ng) % track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Nfloats(ng)/))
# endif
      END IF

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE new_fish
!
#endif
      END MODULE new_fish_mod
