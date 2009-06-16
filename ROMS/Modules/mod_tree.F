      MODULE mod_tree
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
! Because we create and destroy this thing each day, we need an
! efficient way to deallocate everything, hence the need for the gc_list.
!
        USE mod_kinds

        implicit none

        type treenode
          type(treenode), pointer :: left => null()
          type(treenode), pointer :: right => null()
          type(treenode), pointer :: parent => null()
          logical :: red = .FALSE.
          integer :: eggs = 0
          integer :: momfish = 0
          integer :: egg_sum
          real(r8) :: dist
        end type treenode

        type listnode
          type(listnode), pointer :: next => null()
          type(treenode), pointer :: mytreenode => null()
        end type listnode

        type egg_batch
          type(egg_batch), pointer :: next => null()
          type(treenode), pointer :: egg_tree => null()
          type(treenode), pointer :: top => null()
          integer :: egg_sum = 0
        end type egg_batch

        type(treenode), pointer :: tree
        type(listnode), pointer :: gc_list, gc_next
        integer :: Nspawners

        PRIVATE
        PUBLIC :: init, insert, collect, traverse, destroy

        CONTAINS

        SUBROUTINE init
          ALLOCATE(tree)
          ALLOCATE(gc_list)
          gc_list % mytreenode => tree
          gc_next => gc_list
          Nspawners = 0
        END SUBROUTINE init

        SUBROUTINE insert(eggs, dist, ifish)
          integer, intent(in) :: eggs
          integer, intent(in) :: ifish
          real(r8), intent(in) :: dist
          type(treenode), pointer :: cur, p, x, y
          logical :: choice

          ALLOCATE(cur)
          ALLOCATE(gc_next % next)
          gc_next => gc_next % next
          gc_next % mytreenode => cur

          cur % eggs = eggs
          cur % egg_sum = eggs
          cur % dist = dist
          cur % momfish = ifish
          Nspawners = Nspawners + 1

! Empty tree, deposit eggs at the top
          IF (.not. ASSOCIATED(tree % left)) THEN
            tree % left => cur
            cur % parent => tree
            RETURN
          END IF

! Otherwise find somewhere to put these eggs
! New nodes end up at the bottom until a rebalance
          p => tree % left

          DO
            IF (dist <= p % dist) THEN   
              IF (ASSOCIATED(p % left)) THEN
                p => p % left
                CYCLE
              ELSE
                p % left => cur
                cur % parent => p
                EXIT
              END IF
            ELSE
              IF (ASSOCIATED(p % right)) THEN
                p => p % right
                CYCLE
              ELSE
                p % right => cur
                cur % parent => p
                EXIT
              END IF
            END IF
          END DO
! Fix the egg sums above us
          CALL sum_up(cur)
! Balance the thing... red-black for now, until I get smarter about
! balancing eggs.
          cur % red = .true.
          x => cur
          DO WHILE (x % parent % red)
            choice = .false.
            IF (ASSOCIATED(x % parent % parent % left)) THEN
              IF (x % parent % momfish .eq.                             &
     &            x % parent % parent % left % momfish) THEN
                choice = .true.
                IF (ASSOCIATED(x % parent % parent % right)) THEN
                  y => x % parent % parent % right   ! uncle
                  IF (y % red) THEN
                    x % parent % red = .false.
                    y % red = .false.
                    x % parent % parent % red = .true.
                    x => x % parent % parent
                  END IF
                ELSE
                  IF (ASSOCIATED(x % parent % right)) THEN
                    IF (x % momfish .eq.                                &
     &                     x % parent % right % momfish) THEN
!                   IF (x .eq. x % parent % right ) THEN
                      CALL rotate_left(x % parent)
                    END IF
                  END IF
                  x % parent % red = .false.
                  x % parent % parent % red = .true.
                  CALL rotate_right(x % parent % parent)
                END IF
              END IF
            END IF
            IF (ASSOCIATED(x % parent % parent % right) .and. .not.     &
     &                     choice) THEN
              IF (x % parent % momfish .eq.                             &
     &           x % parent % parent % right % momfish) THEN
                IF (ASSOCIATED(x % parent % parent % left)) THEN
                  y => x % parent % parent % left    ! aunt
                  IF (y % red) THEN
                    x % parent % red = .false.
                    y % red = .false.
                    x % parent % parent % red = .true.
                    x => x % parent % parent
                  END IF
                ELSE
                  IF (ASSOCIATED(x % parent % left)) THEN
                    IF (x % momfish .eq.                                &
     &                  x % parent % left % momfish) THEN
                      CALL rotate_right(x % parent)
                    END IF
                  END IF
                  x % parent % red = .false.
                  x % parent % parent % red = .true.
                  CALL rotate_left(x % parent % parent)
                END IF
              END IF
            END IF
          END DO
          tree % left % red = .false.
        END SUBROUTINE insert

! Fix the sums in the ancestor nodes
        SUBROUTINE sum_up(cur)
          type(treenode), pointer :: cur
          type(treenode), pointer :: p
          integer  :: eleft, eright

          p => cur % parent
          DO WHILE (ASSOCIATED(p))
            IF (ASSOCIATED(p % left)) THEN
              eleft = p % left % egg_sum
            ELSE
              eleft = 0
            END IF
            IF (ASSOCIATED(p % right)) THEN
              eright = p % right % egg_sum
            ELSE
              eright = 0
            END IF
            p % egg_sum = eleft + eright + p % eggs
            p => p % parent
          END DO
        RETURN
        END SUBROUTINE sum_up
 
! Fix the sums in me and my parent after a rotation
        SUBROUTINE sum_us(cur)
          type(treenode), pointer :: cur
          type(treenode), pointer :: p
          integer  :: eleft, eright, i

          p => cur
          DO i=1,2
            IF (ASSOCIATED(p % left)) THEN
              eleft = p % left % egg_sum
            ELSE
              eleft = 0
            END IF
            IF (ASSOCIATED(p % right)) THEN
              eright = p % right % egg_sum
            ELSE
              eright = 0
            END IF
            p % egg_sum = eleft + eright + p % eggs
            p => p % parent
          END DO
        RETURN
        END SUBROUTINE sum_us
!
! For the rotating, I'm using C code I found online for red-black trees,
! with reference to Introduction to Algorithms by Cormen, Leiserson,
! Rivest (Chapter 14). It makes right child of x into the parent of x.
!
        SUBROUTINE rotate_left(x)
          type(treenode), pointer :: x, y, z
          integer :: mine, theirs

          y => x % right

          IF (ASSOCIATED(y % left)) THEN
            x % right => y % left
            x % right % parent => x
          ELSE
            NULLIFY(x % right)
          END IF

          z => x % parent
          y % parent => z

          IF (ASSOCIATED(x % parent % left)) THEN
            IF (x % momfish .eq. x % parent % left % momfish) THEN
              x % parent % left => y
            END IF
          ELSE
            x % parent % right => y
          END IF
          y % left => x
          x % parent => y
          CALL sum_us(x)
        END SUBROUTINE rotate_left

        SUBROUTINE rotate_right(x)
          type(treenode), pointer :: x, y, z

          y => x % left

          IF (ASSOCIATED(y % right)) THEN
            x % left => y % right
            x % left % parent => x
          ELSE
            NULLIFY(x % left)
          END IF

          z => x % parent
          y % parent => z

          IF (ASSOCiATED(x % parent % left)) THEN
            IF (x % momfish .eq. x % parent % left % momfish) THEN
              x % parent % left => y
            END IF
          ELSE
            x % parent % right => y
          END IF
          y % right => x
          x % parent => y
          CALL sum_us(x)
        END SUBROUTINE rotate_right

        SUBROUTINE collect(Nsuper, Nfound, eggs, momfish)
          integer, intent(in)  :: Nsuper
          integer, intent(out) :: Nfound
          integer, intent(out) :: eggs(Nsuper)
          integer, intent(out) :: momfish(Nsuper)
          type(egg_batch), pointer :: l, r, p
          type(egg_batch), pointer :: egg_head
!
! I'm going to keep a sorted linked list of my partial batches,
! splitting them until I have enough to fill my superindividuals or I
! run out of mothers with eggs. Actually, I should count up the moms
! and if there aren't too many, just use all the superindividuals.
!
          Nfound = 0
          IF (Nspawners .eq. 0) THEN
            RETURN
          ELSE IF (Nsuper <= Nspawners) THEN
            IF (ASSOCIATED(tree % left)) CALL                           &
     &        copy_traverse(tree % left, Nsuper, Nfound, eggs, momfish)
            IF (Nfound .ne. Nspawners) print *, "Oh, noes!", Nfound,    &
     &            Nspawners
          ELSE
            ALLOCATE(p)
            p % egg_tree => tree % left
          END IF
        END SUBROUTINE collect

        SUBROUTINE traverse
          IF (ASSOCIATED(tree % left)) CALL traverse_print(tree % left)
        END SUBROUTINE traverse

        RECURSIVE SUBROUTINE copy_traverse(p, Nsuper, Nfound, eggs,     &
     &                   momfish)
          type(treenode), pointer :: p
          integer, intent(in)  :: Nsuper
          integer, intent(inout) :: Nfound
          integer, intent(inout) :: eggs(Nsuper)
          integer, intent(inout) :: momfish(Nsuper)

          IF (ASSOCIATED(p % left)) CALL copy_traverse(p % left,        &
     &                Nsuper, Nfound, eggs, momfish)
          Nfound = Nfound + 1
          eggs(Nfound) = p % eggs
          momfish(Nfound) = p % momfish
          IF (ASSOCIATED(p % right)) CALL copy_traverse(p % right,      &
     &                Nsuper, Nfound, eggs, momfish)
        END SUBROUTINE copy_traverse

        RECURSIVE SUBROUTINE traverse_print(p)
          type(treenode), pointer :: p

          IF (ASSOCIATED(p % left)) CALL traverse_print(p % left)
          CALL print_node(p)
          IF (ASSOCIATED(p % right)) CALL traverse_print(p % right)

        END SUBROUTINE traverse_print

        SUBROUTINE print_node(cur)
          type(treenode), pointer :: cur
          print *, "Node: ", cur % dist, cur % eggs, cur % egg_sum
        END SUBROUTINE print_node

        SUBROUTINE destroy
! Loop through, cleaning up both the treenodes and the listnodes.
          type(listnode), pointer :: cur, last

          cur => gc_list
          DO WHILE (ASSOCIATED(cur % next))
            DEALLOCATE(cur % mytreenode)
            last => cur
            cur => cur % next
            DEALLOCATE(last)
          END DO
          DEALLOCATE(cur % mytreenode)
          DEALLOCATE(cur)
        END SUBROUTINE destroy

      END MODULE mod_tree
