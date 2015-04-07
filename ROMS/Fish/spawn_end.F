#include "cppdefs.h"
      MODULE spawn_end
#if defined NEMURO_SAN
!
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
        USE mod_kinds
        USE mod_types
        USE mod_param
        USE mod_biology

        implicit none

        type fish_node
          type(fish_node), pointer :: next => null()
          integer :: fish
          real(r8) :: fish_worth
        end type fish_node

!...........................................................................
! Global nil node for leaves and parent of root.
!...........................................................................
        type(fish_node), target, save :: nil

        type(fish_node), pointer :: fish_list

        PRIVATE :: push_fish, pop_fish, clone_fish, clean_list,         &
     &             isLast, nil
        PUBLIC  :: split_biggest

        CONTAINS

        SUBROUTINE split_biggest(ng, isp, Navail)
          USE mod_fish
          integer, intent(in) :: ng
          integer, intent(in) :: isp
          integer, intent(inout) :: Navail

          integer :: ifish, nfish, nf
          type(fish_node), pointer :: this_fish
          type(fish_node), pointer :: half_fish
          real(r8) :: worth

!  Store all the first-year fish of this species
          ALLOCATE(fish_list)
          fish_list % next => nil
          nf = 0
          DO ifish = FISHES(ng)%age_base(isp)+1,                        &
     &              FISHES(ng)%age_base(isp)+Nfishperyear(ng)
! use current worth to split individuals
!            worth = FISHES(ng)%bioenergy(ifworth,ifish)
! use initial worth to split individuals
            worth = FISHES(ng)%bioenergy(ifiniwth,ifish)
!           print *, 'fish x has worth y', ifish, worth
            IF (worth .gt. 0.1_r8) THEN
              ALLOCATE(this_fish)
              this_fish % fish = ifish
              this_fish % fish_worth = worth
              CALL push_fish(this_fish)
              nf = nf + 1
            ELSE
              EXIT
            END IF
          END DO

!  Find the largest and split them in two
          DO WHILE (Navail .gt. 0)
            this_fish => pop_fish()
            this_fish % fish_worth = 0.5 * this_fish % fish_worth
            ifish = this_fish % fish
            nfish = FISHES(ng) % next_free(isp)
            FISHES(ng) % next_free(isp) =                               &
     &                  FISHES(ng) % next_free(isp) + 1
            ALLOCATE(half_fish)
            half_fish % fish = nfish
            half_fish % fish_worth = this_fish % fish_worth
            CALL clone_fish(ng, ifish, nfish)
            CALL push_fish(this_fish)
            CALL push_fish(half_fish)
            Navail = Navail - 1
          END DO
          CALL clean_list
        END SUBROUTINE split_biggest

        FUNCTION pop_fish()
          type(fish_node), pointer :: pop_fish
          type(fish_node), pointer :: cur, prev

          prev => fish_list
!  Returning the first box.
          cur => fish_list % next
          prev % next => cur % next
          pop_fish => cur
          RETURN
        END FUNCTION pop_fish

        SUBROUTINE push_fish(box)
          type(fish_node), pointer :: box
          type(fish_node), pointer :: cur, next
          real(r8) :: cur_fish, next_fish, my_fish

          cur => fish_list
          next => cur % next
          my_fish = box % fish_worth
          cur_fish = 1.e35
          IF (.not. isLast(next)) THEN
            next_fish = next % fish_worth
          ELSE
            next_fish = 0
          END IF
          DO
            IF (my_fish < cur_fish .and. my_fish >= next_fish) THEN
              box % next => next
              cur % next => box
              RETURN
            END IF
            cur => next
            next => cur % next
            cur_fish = cur % fish_worth
            IF (.not. isLast(next)) THEN
              next_fish = next % fish_worth
            ELSE
              next_fish = 0
            END IF
          END DO
        END SUBROUTINE push_fish

        SUBROUTINE clone_fish(ng, fish1, fish2)
        USE mod_param 
        USE mod_fish
        USE mod_scalars
        USE nrutil
        USE ran_state, ONLY: ran_seed

        integer, intent(in) :: ng, fish1, fish2

        integer :: i, j, Ir, Jr
        real(r8) :: xnudg, ynudg, xfish, yfish

!  Copy location and other features of first fish

          print*, 'CLONE', fish1, fish2,                                &
     &                     FISHES(ng)%bioenergy(ifiniwth,fish1)

          FISHES(ng)%bioenergy(ifworth,fish1)=0.5_r8*                   &
     &                           FISHES(ng)%bioenergy(ifworth,fish1)
          FISHES(ng)%bioenergy(ifiniwth,fish1)=0.5_r8*                   &
     &                           FISHES(ng)%bioenergy(ifiniwth,fish1)
          FISHES(ng)%egg_num(fish1)=0.5_r8*FISHES(ng)%egg_num(fish1)
          FISHES(ng)%ysac_num(fish1)=0.5_r8*FISHES(ng)%ysac_num(fish1)
          FISHES(ng)%larv_num(fish1)=0.5_r8*FISHES(ng)%larv_num(fish1)
          FISHES(ng)%juv_num(fish1)=0.5_r8*FISHES(ng)%juv_num(fish1)

          FISHES(ng) % bounded(fish2) = .TRUE.

          DO j=0,NFT
            DO i=1,NFV(ng)
              FISHES(ng) % track(i,j,fish2) =                           &
     &              FISHES(ng) % track(i,j,fish1)
            END DO
          END DO
          DO i=1,NFishV(ng)
            FISHES(ng) % bioenergy(i,fish2) =                           &
     &              FISHES(ng) % bioenergy(i,fish1)
          END DO
          FISHES(ng)%egg_dur(fish2) = FISHES(ng)%egg_dur(fish1)
          FISHES(ng)%egg_num(fish2) = FISHES(ng)%egg_num(fish1)
          FISHES(ng)%ysac_dur(fish2)= FISHES(ng)%ysac_dur(fish1)
          FISHES(ng)%ysac_num(fish2)= FISHES(ng)%ysac_num(fish1)
          FISHES(ng)%larv_dur(fish2)= FISHES(ng)%larv_dur(fish1)
          FISHES(ng)%larv_num(fish2)= FISHES(ng)%larv_num(fish1)
          FISHES(ng)%juv_dur(fish2) = FISHES(ng)%juv_dur(fish1)
          FISHES(ng)%juv_num(fish2) = FISHES(ng)%juv_num(fish1)
          FISHES(ng)%lifestage(fish2) = FISHES(ng)%lifestage(fish1)
          FISHES(ng)%swimtype(1,fish2) = FISHES(ng)%swimtype(1,fish1)
          FISHES(ng)%swimtype(2,fish2) = FISHES(ng)%swimtype(2,fish1)
          FISHES(ng)%species(fish2) = FISHES(ng)%species(fish1)
          FISHES(ng)%alive(fish2) = FISHES(ng)%alive(fish1)

!  Randomize the location just a little
          CALL ran1 (xnudg)
          CALL ran1 (ynudg)
          xfish = FISHES(ng)%track(ixgrd,0,fish2)
          yfish = FISHES(ng)%track(iygrd,0,fish2)
          DO j=0,NFT
            FISHES(ng)%track(ixgrd,j,fish2) = xfish + (xnudg-0.5_r8)
            FISHES(ng)%track(iygrd,j,fish2) = yfish + (ynudg-0.5_r8)
          END DO
          FISHES(ng)%Tinfo(ixgrd,fish2)=FISHES(ng)%Tinfo(ixgrd,fish1)
          FISHES(ng)%Tinfo(iygrd,fish2)=FISHES(ng)%Tinfo(iygrd,fish1)
          FISHES(ng)%Tinfo(izgrd,fish2)=FISHES(ng)%Tinfo(izgrd,fish1)
          Ir=NINT(FISHES(ng)%track(ixgrd,0,fish2))
          Jr=NINT(FISHES(ng)%track(iygrd,0,fish2))
          FISHES(ng)%cellid(fish2)=Ir+(Lm(ng)+2)*Jr

        END SUBROUTINE clone_fish

!...........................................................................
        FUNCTION isLast(x) result(b)
!..........................................................................
! Check if node x is the last one
!..........................................................................
          type (fish_node), pointer :: x
          logical :: b

          b = ASSOCIATED(x, nil)
        END FUNCTION isLast

        SUBROUTINE clean_list
          type(fish_node), pointer :: cur, next

          cur => fish_list % next
          DEALLOCATE(fish_list)
          DO WHILE (.not. isLast(cur))
            next => cur % next
!           print *, 'deallocating', cur % fish_worth, cur % fish
            DEALLOCATE(cur)
            cur => next
          END DO
        END SUBROUTINE clean_list
#endif
      END MODULE spawn_end
