      MODULE mod_trees
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
        implicit none

        type treenode
          type(treenode), pointer :: left => null()
          type(treenode), pointer :: right => null()
          type(treenode), pointer :: parent => null()
          logical :: red = .FALSE.
          logical :: leaf
          integer :: eggs
          real*8  :: dist
        end type treenode

        type(treenode), pointer :: tree

        PRIVATE
        PUBLIC :: init, insert, traverse, destroy

        CONTAINS

        SUBROUTINE init
          ALLOCATE(tree)
        END SUBROUTINE init

        SUBROUTINE insert(eggs, dist)
          integer, intent(in) :: eggs
          real*8, intent(in)  :: dist
          type(treenode), pointer :: cur, p, new

          ALLOCATE(cur)
          cur % eggs = eggs
          cur % dist = dist
          cur % leaf = .TRUE.

! Empty tree, go left and deposit eggs
          IF (.not. ASSOCIATED(tree % left)) THEN
            tree % left => cur
            tree % dist = dist
            tree % eggs = eggs
            tree % leaf = .FALSE.
            cur % parent => tree
            RETURN
          END IF

! Otherwise find somewhere to put these eggs
          p => tree

          DO
            IF (dist < p % dist) THEN   
              IF (ASSOCIATED(p % left)) THEN
                p => p % left
                CYCLE
              ELSE IF (p % leaf) THEN
                ALLOCATE(new)
                new % leaf = .true.
                new % eggs = p % eggs
                new % dist = p % dist
                new % parent => p
                p % right => new
                p % leaf = .false.
                p % left => cur
                cur % parent => p
                CALL sum_up(cur)
              ELSE
                p % left => cur
                cur % parent => p
                CALL sum_up(cur)
                EXIT
              END IF
            ELSE
              IF (ASSOCIATED(p % right)) THEN
                p => p % right
                CYCLE
              ELSE
                p % right => cur
                p % dist = dist
                p % eggs = eggs
                cur % parent => p
                CALL sum_up(p)
                EXIT
              END IF
            END IF
          END DO
        END SUBROUTINE insert

! Insertion messes up the sums in the ancestor nodes
        SUBROUTINE sum_up(cur)
          type(treenode), pointer :: cur
          type(treenode), pointer :: p
          integer :: eleft, eright
          real*8  :: dleft, dright

          p => cur % parent
          DO WHILE (ASSOCIATED(p))
            IF (ASSOCIATED(p % left)) THEN
              eleft = p % left % eggs
              dleft = p % left % dist
            ELSE
              eleft = 0
            END IF
            IF (ASSOCIATED(p % right)) THEN
              eright = p % right % eggs
              dright = p % right % dist
            ELSE
              eright = 0
            END IF
            p % eggs = eleft + eright
            p % dist = (eleft*dleft + eright*dright) / (eleft + eright)
            p => p % parent
          END DO
          
        END SUBROUTINE sum_up

      END MODULE mod_trees
