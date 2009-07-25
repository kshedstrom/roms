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
! What I'm trying here is to build a balanced binary tree of the eggs.
! Each node represents one spawning female. It needs to know how many
! eggs, so as to group them into the resulting superindividuals allowed
! for that day. For the grouping, it needs some "distance", be it
! longshore distance or some other idealized metric.
!
! The distance is "spawn_dist" on the model grid. A spawning female can
! look this up based on its position in i,j space.
!
! Idea for the future - instead of one treetop, have one treetop per
! species per region. Regions could be divided by domains in spawn_dist,
! i.e., have the 1000's region, the 2000's region, and so on. If one
! species had three regions, it would require at least three
! superindividuals to fill per day: load the roe_list accordingly.
!
! Because we create and destroy this thing each day, we need an
! efficient way to deallocate everything; hence the need for the gc_list.
!
        USE mod_kinds
        USE mod_types

        implicit none

        type egg_box
	  type(egg_box), pointer :: next => null()
          integer :: i1, i2, j1, j2
          integer :: egg_sum
        end type egg_box

        type(egg_box), pointer :: eggbox_head

        PRIVATE
        PUBLIC :: egglist_init, egglist_collect, egglist_destroy

        CONTAINS

! This insert is visible from outside the module
        SUBROUTINE egglist_init
          ALLOCATE(eggbox_head)
          Nspawners = 0
        END SUBROUTINE egglist_init

        SUBROUTINE egglist_collect(Nsuper, Nfound, eggs, momfish)
          integer, intent(in)  :: Nsuper
          integer, intent(out) :: Nfound
          real(r8), intent(out) :: eggs(Nsuper)
          integer, intent(out) :: momfish(Nsuper)
          type(treenode), pointer :: l_tree, r_tree, t_tree
          logical :: single
!
! First, count up the moms and if there aren't too many, just use all
! the superindividuals.
!
          Nfound = 0
          IF (Nspawners .eq. 0) THEN
            RETURN
          ELSE IF (Nsuper >= Nspawners) THEN
            IF (ASSOCIATED(treetop % left)) CALL copy_traverse(         &
     &               treetop % left, Nsuper, Nfound, eggs, momfish)
            IF (Nfound .ne. Nspawners) print *, "Oh, noes!", Nfound,    &
     &            Nspawners
          ELSE
!
! Now we have to group the eggs into the available superindividuals.
! I'm going to keep a sorted linked list of my partial batches,
! splitting them until I have enough to fill my superindividuals.
!
            ALLOCATE(roe_list % next)
            ALLOCATE(gc_next % next)
            gc_next => gc_next % next
            gc_next % myeggs => roe_list % next
            roe_list % next % egg_tree => treetop
            roe_list % next % egg_sum = treetop % left % egg_sum
            Nfound = 1

            DO WHILE (Nfound .lt. Nsuper)
! Pop the top tree from the stack and split it.
              single = .false.
              t_tree => pop_tree(.true.)
              IF (.not. ASSOCIATED(t_tree)) THEN
                print *, "I'm in trouble again..."
              END IF

! Need to make some new treenodes: these become the "treetop" parts of
! the left and right subtrees, above the actual egg-carrying nodes.
              ALLOCATE(l_tree)
              ALLOCATE(r_tree)
              ALLOCATE(gc_next % next)
              gc_next => gc_next % next
              gc_next % mytreenode => l_tree
              ALLOCATE(gc_next % next)
              gc_next => gc_next % next
              gc_next % mytreenode => r_tree

              IF (ASSOCIATED (t_tree % left)) THEN
                l_tree % left => t_tree % left
                nullify(t_tree % left)
                l_tree % egg_sum = l_tree % left % egg_sum
              ELSE
                l_tree % left => t_tree
                single = .true.
                t_tree % egg_sum = t_tree % eggs
                l_tree % egg_sum = t_tree % eggs
              END IF
              IF (ASSOCIATED (t_tree % right)) THEN
                r_tree % left => t_tree % right
                nullify(t_tree % right)
                r_tree % egg_sum = r_tree % left % egg_sum
              ELSE
                r_tree % left => t_tree
                single = .true.
                t_tree % egg_sum = t_tree % eggs
                r_tree % egg_sum = t_tree % eggs
              END IF
              l_tree % left % parent => l_tree
              l_tree % left % red = .false.
              r_tree % left % parent => r_tree
              r_tree % left % red = .false.
! If we've both subtrees, give the top's eggs to the smaller
! (by egg count) tree, making a new treenode (because that's what
! insert_priv does).
              IF (.not. single) THEN
                IF (l_tree % left % egg_sum <                           &
     &                     r_tree % left % egg_sum) THEN
                  CALL insert_priv(t_tree % eggs, t_tree % dist,        &
     &                 t_tree % momfish, l_tree)
                ELSE
                  CALL insert_priv(t_tree % eggs, t_tree % dist,        &
     &                 t_tree % momfish, r_tree)
                END IF
              END IF
              CALL push_trees(l_tree, r_tree)
              Nfound = Nfound + 1
            END DO
            Nfound = 0
            DO WHILE (Nfound .lt. Nsuper)
              t_tree => pop_tree(.false.)
              IF (ASSOCIATED(t_tree)) THEN
                Nfound = Nfound + 1
                eggs(Nfound) = t_tree % egg_sum 
                momfish(Nfound) = t_tree % momfish  ! one of the bunch
              ELSE
                RETURN
              END IF
            END DO
          END IF
        END SUBROUTINE egglist_collect

        SUBROUTINE egglist_destroy
! Loop through, cleaning up both the treenodes and the listnodes.
          type(listnode), pointer :: cur, last

          cur => gc_list
          DO WHILE (ASSOCIATED(cur % next))
            IF (ASSOCIATED(cur % mytreenode))                           &
     &          DEALLOCATE(cur % mytreenode)
            IF (ASSOCIATED(cur % myeggs))                               &
     &          DEALLOCATE(cur % myeggs)
            last => cur
            cur => cur % next
            DEALLOCATE(last)
          END DO
          DEALLOCATE(cur % mytreenode)
          DEALLOCATE(cur)
        END SUBROUTINE egglist_destroy

#endif
      END MODULE mod_egglist
