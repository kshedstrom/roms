#include "cppdefs.h"
      MODULE mod_egglist
#if defined NEMURO_SAN && defined EGGS_BISECTION
!
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!  Set up tree structure and functions.                                !
!=======================================================================
!
! This is Chris Edward's "bisection" idea. Start by knowing the locations
! of all the eggs. The routine eggs_by_cell finds the eggs for each cell,
! but only for that tile. We start by collecting them into one big array
! on the master.
!
! Because we create and destroy this thing each day, we need an
! efficient way to deallocate everything; hence the need for the gc_list.
!
        USE mod_kinds
        USE mod_types
        USE mod_param

        implicit none

        type egg_box
          type(egg_box), pointer :: next => null()
          integer :: i1, i2, j1, j2
          integer :: egg_sum
          logical :: last_lr_split
        end type egg_box

!...........................................................................
! Global nil node for leaves and parent of root.
!...........................................................................
        type(egg_box), target, save :: nil

        real(r8), allocatable :: eggs_all(:,:,:)
        type(egg_box), pointer :: roe_list(:)

        PRIVATE :: eggs_all, push_box, pop_box, isLast, find_fish, nil
        PUBLIC  :: eggs_share, egglist_init, egglist_collect,           &
     &             egglist_destroy, egglist_split

        CONTAINS

        SUBROUTINE eggs_share(ng, model, LBi, UBi, LBj, UBj, egg_count)
          USE mod_parallel
# ifdef DISTRIBUTE
          USE distribute_mod, ONLY : mp_bcasti, mp_gather2d
# endif
          integer, intent(in) :: ng, model
          integer, intent(in) :: LBi, UBi, LBj, UBj 
# ifdef ASSUMED_SHAPE
          real(r8), intent(in) :: egg_count(LBi:,LBj:,:)
# else
          real(r8), intent(in) ::                                       &
     &              egg_count(LBi:UBi,LBj:UBj,Nspecies(ng))
# endif
          integer :: isp, i, j
# ifdef DISTRIBUTE
          real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*isp) :: eggs_line
# endif
          
! First, collect every tile's eggs into global array
# ifdef DISTRIBUTE
          eggs_line = 0.0_r8
# endif
          IF (.not. ALLOCATED(eggs_all)) THEN
            allocate(eggs_all(0:(Lm(ng)+1),0:(Mm(ng)+1),isp))
          END IF
# ifdef DISTRIBUTE
          CALL mp_gather3d (ng, model, LBi, UBi, LBj, UBj,              &
     &                  1, Nspecies(ng), 1, r3dvar, 1.0_r8,             &
#  ifdef MASKING
     &                  Amask,                                          &
#  endif
     &                  egg_count, Npts, eggs_line, SetFillVal)
        eggs_all = reshape(eggs_line,(/ Lm(ng)+2, Mm(ng)+2, Nspecies /))
# else
        DO isp=1,Nspecies(ng)
          DO j=0,Mm(ng)+1
            DO i=0,Mm(ng)+1
              eggs_all(i,j,isp)=egg_count(i,j,isp)
            END DO
          END DO
        END DO
# endif
        END SUBROUTINE eggs_share

        SUBROUTINE egglist_init(ng)
          integer, intent(in) :: ng

          integer :: isp
          type(egg_box), pointer :: box

          ALLOCATE(roe_list(Nspecies(ng)))

          DO isp=1,Nspecies(ng)
            ALLOCATE(box)
            roe_list(isp) % next => box
            box % next => nil
            box % i1 = 0
            box % i2 = Lm(ng)+1
            box % j1 = 0
            box % j2 = Mm(ng)+1
            IF (Lm(ng) > Mm(ng)) THEN
              box % last_lr_split = .false.
            ELSE
              box % last_lr_split = .true.
            END IF
            box % egg_sum = SUM(eggs_all(:,:,isp))
          END DO
        END SUBROUTINE egglist_init

        SUBROUTINE egglist_split(isp, Nsuper)
          integer, intent(in)  :: isp
          integer, intent(in)  :: Nsuper
          type(egg_box), pointer :: bigbox, ur_box, ll_box
          integer :: Nboxes
!
! First, if there's only one superindividual, put all the eggs in one
! basket.
!
          Nboxes = 1
          bigbox => roe_list(isp) % next
          IF (Nsuper .eq. 0) THEN
            IF (bigbox % egg_sum .ne. 0) THEN
              print *, 'TROUBLE: eggs without a home ', isp,            &
     &                bigbox % egg_sum
            END IF
            RETURN
          ELSE IF (Nsuper .eq. 1) THEN
            RETURN
          ELSE
!
! Now we have to split the eggs into the available superindividuals.
! I'm going to keep a sorted linked list of my boxes, splitting
! them until I have enough to fill my superindividuals.
!
            DO WHILE (Nboxes .lt. Nsuper)
! Pop the top box from the stack and split it.
              bigbox => pop_box(isp)
              IF (.not. ASSOCIATED(bigbox)) THEN
                print *, "I'm in trouble again..."
              END IF

! Need to make some new boxes.
              ALLOCATE(ll_box)
              ALLOCATE(ur_box)

              IF (bigbox % last_lr_split) THEN
!  split north-south
                ll_box % i1 = bigbox % i1
                ll_box % i2 = bigbox % i2
                ll_box % j1 = bigbox % j1
                ll_box % j2 = (bigbox % j1 + bigbox % j2)/2
                ur_box % i1 = bigbox % i1
                ur_box % i2 = bigbox % i2
                ur_box % j1 = (bigbox % j1 + bigbox % j2)/2 + 1
                ur_box % j2 = bigbox % j2
                ll_box % last_lr_split = .false.
                ur_box % last_lr_split = .false.
              ELSE
!  split left-right
                ll_box % i1 = bigbox % i1
                ll_box % i2 = (bigbox % i1 + bigbox % i2)/2
                ll_box % j1 = bigbox % j1
                ll_box % j2 = bigbox % j2
                ur_box % i1 = (bigbox % i1 + bigbox % i2)/2 + 1
                ur_box % i2 = bigbox % i2
                ur_box % j1 = bigbox % j1
                ur_box % j2 = bigbox % j2
                ll_box % last_lr_split = .false.
                ur_box % last_lr_split = .false.
              END IF
              Nboxes = Nboxes - 1
              ll_box % egg_sum = SUM(eggs_all(ll_box % i1:ll_box % i2, &
     &                           ll_box % j1:ll_box % j2,isp))
              ur_box % egg_sum = SUM(eggs_all(ur_box % i1:ur_box % i2, &
     &                           ur_box % j1:ur_box % j2,isp))
              IF (ll_box % egg_sum > 0) THEN
                CALL push_box(isp, ll_box)
                Nboxes = Nboxes + 1
              ELSE
                deallocate(ll_box)
              END IF
              IF (ur_box % egg_sum > 0) THEN
                CALL push_box(isp, ur_box)
                Nboxes = Nboxes + 1
              ELSE
                deallocate(ur_box)
              END IF
              deallocate(bigbox)
            END DO
          END IF
        END SUBROUTINE egglist_split

        SUBROUTINE egglist_collect(isp, Nsuper, Nfound, eggs,           &
     &             ifish, jfish)
          integer, intent(in)  :: isp
          integer, intent(in)  :: Nsuper
          integer, intent(out) :: Nfound
          real(r8), intent(out) :: eggs(Nsuper)
          integer, intent(out) :: ifish(Nsuper)
          integer, intent(out) :: jfish(Nsuper)
          type(egg_box), pointer :: box
          logical :: single
!
          Nfound = 0
          DO
            box => pop_box(isp)
            IF (.not. isLast(box)) THEN
              Nfound = Nfound + 1
              eggs(Nfound) = box % egg_sum 
              CALL find_fish(isp, box, ifish(Nfound), jfish(Nfound))
            ELSE
              RETURN
            END IF
          END DO
        END SUBROUTINE egglist_collect

        FUNCTION pop_box(isp)
          integer, intent(in) :: isp
          type(egg_box), pointer :: pop_box
          type(egg_box), pointer :: cur, prev

          prev => roe_list(isp)
!  Returning the first box.
          cur => roe_list(isp) % next
          prev % next => cur % next
          pop_box => cur
          RETURN
        END FUNCTION pop_box

        SUBROUTINE push_box(isp, box)
          integer, intent(in) :: isp
          type(egg_box), pointer :: box
          type(egg_box), pointer :: cur, next
          real(r8) :: cur_eggs, next_eggs, my_eggs

          cur => roe_list(isp)
          next => cur % next
          my_eggs = box % egg_sum
          cur_eggs = 1.e35
          IF (.not. isLast(next)) THEN
            next_eggs = next % egg_sum
          ELSE
            next_eggs = 0
          END IF
          DO
            IF (my_eggs < cur_eggs .and. my_eggs >= next_eggs) THEN
              box % next => next
              cur % next => box
              RETURN
            END IF
            cur => next
            next => cur % next
            cur_eggs = cur % egg_sum
            IF (.not. isLast(next)) THEN
              next_eggs = next % egg_sum
            ELSE
              next_eggs = 0
            END IF
          END DO
        END SUBROUTINE push_box

!...........................................................................
        FUNCTION isLast(x) result(b)
!..........................................................................
! Check if node x is the last one
!..........................................................................
          type (egg_box), pointer :: x
          logical :: b

          b = ASSOCIATED(x, nil)
        END FUNCTION isLast

        SUBROUTINE egglist_destroy(isp)
        integer, intent(in) :: isp

          DEALLOCATE(roe_list)
        END SUBROUTINE egglist_destroy

        SUBROUTINE find_fish(isp, box, ifish, jfish)
          type(egg_box), pointer :: box
          integer, intent(in)  :: isp
          integer, intent(out) :: ifish, jfish
          real(r8) :: egg_max
          integer :: i, j

          ifish = box%i1
          jfish = box%j1
          egg_max = eggs_all(box%i1, box%j1, isp)
          DO j=box%j1, box%j2
            DO i=box%i1, box%i2
              IF (eggs_all(i,j,isp) > egg_max) THEN
                ifish = i
                jfish = j
                egg_max = eggs_all(i,j,isp)
              END IF
            END DO
          END DO
        END SUBROUTINE find_fish

#endif
      END MODULE mod_egglist
