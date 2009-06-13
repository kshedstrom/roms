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
! Each leaf node represents one spawning female. It needs to know how
! many eggs, so as to group them into the resulting eggs objects allowed
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
          integer :: egg_sum
          real(r8) :: dist
        end type treenode

        type listnode
          type(listnode), pointer :: next => null()
          type(treenode), pointer :: mytreenode => null()
        end type listnode

        type(treenode), pointer :: tree
        type(listnode), pointer :: gc_list, gc_next

        PRIVATE
        PUBLIC :: init, insert, traverse, destroy

        CONTAINS

        SUBROUTINE init
          ALLOCATE(tree)
          ALLOCATE(gc_list)
          gc_list % mytreenode => tree
          gc_next => gc_list
        END SUBROUTINE init

        SUBROUTINE insert(eggs, dist)
          integer, intent(in) :: eggs
          real(r8), intent(in) :: dist
          type(treenode), pointer :: cur, p, new

! Empty tree, deposit eggs here
          IF (tree % eggs .eq. 0) THEN
            tree % dist = dist
            tree % eggs = eggs
            tree % egg_sum = eggs
            RETURN
          END IF

! New nodes end up at the bottom until a rebalance
          ALLOCATE(cur)
          ALLOCATE(gc_next % next)
          gc_next => gc_next % next
          gc_next % mytreenode => cur

          cur % eggs = eggs
          cur % egg_sum = eggs
          cur % dist = dist

! Otherwise find somewhere to put these eggs
          p => tree

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
! Fix the parent sums
          CALL sum_up(cur)
! Balance the thing...          
          IF (cur % egg_sum > p % egg_sum) THEN
          END IF
        END SUBROUTINE insert

! Insertion messes up the sums in the ancestor nodes
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

        SUBROUTINE traverse
          CALL traverse_priv(tree)
        END SUBROUTINE traverse

        RECURSIVE SUBROUTINE traverse_priv(p)
          type(treenode), pointer :: p

          IF (ASSOCIATED(p % left)) CALL traverse_priv(p % left)
          CALL print_node(p)
          IF (ASSOCIATED(p % right)) CALL traverse_priv(p % right)

        END SUBROUTINE traverse_priv

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
