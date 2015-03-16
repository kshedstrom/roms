      MODULE mod_tree
!
!================================================== Kate Hedstrom ======
!    Fixes and improvements by David Car (david.car7@gmail.com)        !
!=======================================================================
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
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

        implicit none

        type treenode
          type(treenode), pointer :: left => null()
          type(treenode), pointer :: right => null()
          type(treenode), pointer :: parent => null()
          logical :: red = .FALSE.
          real(r8) :: eggs = 0
          integer :: momfish = 0
          real(r8) :: egg_sum
          real(r8) :: dist
        end type treenode

        type RedBlackTree
           type(treenode), pointer :: root => null()
           integer :: nitems = 0
        end type RedBlackTree

!...........................................................................
! Global nil node for leaves and parent of root.
!...........................................................................
        type(treenode), target, save :: nil

        type listnode
          type(listnode), pointer :: next => null()
          type(treenode), pointer :: mytreenode => null()
          type(egg_batch), pointer :: myeggs => null()
        end type listnode

        type egg_batch
          type(egg_batch), pointer :: next => null()
          type(treenode), pointer :: egg_tree => null()
          real(r8) :: egg_sum = 0
        end type egg_batch

        type(listnode), pointer :: gc_list, gc_next
        type(egg_batch), pointer :: roe_list

        PRIVATE :: insert_priv, sum_up, sum_us, rotate_left,            &
     &            rotate_right, isLeaf, isLeft, isRight, pop_tree,      &
     &            push_trees, push_one_tree, copy_traverse,             &
     &            traverse_print, print_node
        PUBLIC :: tree_init, tree_insert, tree_collect, tree_traverse,  &
     &            tree_destroy

        CONTAINS

        SUBROUTINE tree_init(this)
          type (RedBlackTree) :: this
          ALLOCATE(this % root)
          ALLOCATE(gc_list)
          ALLOCATE(roe_list)
          this % root % left => nil
          this % root % right => nil
          this % root % parent => nil
          gc_list % mytreenode => this % root
          gc_list % myeggs => roe_list
          gc_next => gc_list
        END SUBROUTINE tree_init

! This insert is visible from outside the module
        SUBROUTINE tree_insert(this, eggs, dist, ifish)
          type (RedBlackTree) :: this
          real(r8), intent(in) :: eggs
          integer, intent(in) :: ifish
          real(r8), intent(in) :: dist

          IF (.not. ASSOCIATED(this % root)) RETURN
          this % nitems = this % nitems + 1
          CALL insert_priv(this % root, eggs, dist, ifish)

        END SUBROUTINE tree_insert

        SUBROUTINE insert_priv(top, eggs, dist, ifish)
          type (treenode), target :: top
          type (treenode), pointer :: treetop
          real(r8), intent(in) :: eggs
          integer, intent(in) :: ifish
          real(r8), intent(in) :: dist
          type(treenode), pointer :: cur, p, x, y, z

          treetop => top

          ALLOCATE(cur)
          ALLOCATE(gc_next % next)
          gc_next => gc_next % next
          gc_next % mytreenode => cur

          cur % eggs = eggs
          cur % egg_sum = eggs
          cur % dist = dist
          cur % momfish = ifish
          cur % left => nil
          cur % right => nil
          cur % parent => nil

! Empty tree, deposit eggs at the top
         IF (isLeaf(treetop % left)) THEN
           treetop % left => cur
           RETURN
         ENDIF

! Otherwise find somewhere to put these eggs
! New nodes end up at the bottom until a rebalance
          p => treetop % left

          DO
            IF (dist <= p % dist) THEN
              IF (.not. isLeaf(p % left)) THEN
                p => p % left
                CYCLE
              ELSE
                p % left => cur
                cur % parent => p
                EXIT
              END IF
            ELSE
              IF (.not. isLeaf(p % right)) THEN
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
! The balancing isn't working, so let's just return here.
!         RETURN
! Balance the thing... red-black for now.
          cur % red = .true.
          x => cur
          p => null()
          DO WHILE (x % parent % red)
            IF (ASSOCIATED(x % parent % parent % left, x % parent)) THEN
              y => x % parent % parent % right   ! uncle
              IF (y % red) THEN
                x % parent % red = .false.
                y % red = .false.
                x % parent % parent % red = .true.
                x => x % parent % parent
              ELSE
                IF (ASSOCIATED(x, x % parent % right)) THEN
                  x => x % parent
                  CALL rotate_left(treetop, x)
                END IF
                x % parent % red = .false.
                x % parent % parent % red = .true.
                CALL rotate_right(treetop, x % parent % parent)
              END IF
            ELSE         ! parent is grandparent's right child
              y => x % parent % parent % left    ! aunt
              IF (y % red) THEN
                x % parent % red = .false.
                y % red = .false.
                x % parent % parent % red = .true.
                x => x % parent % parent
              ELSE
                IF (ASSOCIATED(x, x % parent % left)) THEN
                  x => x % parent
                  CALL rotate_right(treetop, x)
                END IF
                x % parent % red = .false.
                x % parent % parent % red = .true.
                CALL rotate_left(treetop, x % parent % parent)
              END IF
            END IF
          END DO
          treetop % left % red = .false.
        END SUBROUTINE insert_priv

! Fix the sums in the ancestor nodes
        SUBROUTINE sum_up(cur)
          type(treenode), pointer :: cur
          type(treenode), pointer :: p
          real(r8)  :: eleft, eright

          p => cur % parent
          DO WHILE (.not. ASSOCIATED(p, nil))
            IF (.not. isLeaf(p % left)) THEN
              eleft = p % left % egg_sum
            ELSE
              eleft = 0
            END IF
            IF (.not. isLeaf(p % right)) THEN
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
          real(r8)  :: eleft, eright
          integer  :: i

          p => cur
          DO i=1,2
            IF (.not. isLeaf(p % left)) THEN
              eleft = p % left % egg_sum
            ELSE
              eleft = 0
            END IF
            IF (.not. isLeaf(p % right)) THEN
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
! Rivest (Chapter 14). It makes right child of x_ into the parent of x_.
!
        SUBROUTINE rotate_left(top, x_)
          type(treenode), target, intent(inout) :: top
          type(treenode), target, intent(inout) :: x_
          type(treenode), pointer :: treetop
          type(treenode), pointer :: x => null()
          type(treenode), pointer :: y => null()

          treetop => top
          x => x_
          y => x % right

          x % right => y % left

          IF (.not. isLeaf(y % left)) THEN
             y % left % parent => x
          END IF

          y % parent => x % parent

          IF (ASSOCIATED(x % parent, nil)) THEN
             treetop % left => y
          ELSE
            IF (ASSOCIATED(x, x % parent % left)) THEN
              x % parent % left => y
            ELSE
              x % parent % right => y
            END IF
          END IF
          y % left => x
          x % parent => y
          CALL sum_us(x)
        END SUBROUTINE rotate_left

        SUBROUTINE rotate_right(top, x_)
          type(treenode), target, intent(inout) :: top
          type(treenode), target, intent(inout) :: x_
          type(treenode), pointer :: treetop
          type(treenode), pointer :: x => null()
          type(treenode), pointer :: y => null()

          treetop => top
          x => x_
          y => x % left

          x % left => y % right

          IF (.not. isLeaf(y % right)) THEN
            y % right % parent => x
          END IF

          y % parent => x % parent

          IF (ASSOCIATED(x % parent, nil)) THEN
            treetop % left => y
          ELSE
            IF (ASSOCIATED(x, x % parent % right)) THEN
              x % parent % right => y
            ELSE
              x % parent % left => y
            END IF
          END IF
          y % right => x
          x % parent => y
          CALL sum_us(x)
        END SUBROUTINE rotate_right

!...........................................................................
! None of these are used at the moment
!...........................................................................

!...........................................................................
        FUNCTION isLeft(x, y) result(b)
!..........................................................................
! Check if node y is left child of x
!..........................................................................
          type (treenode), pointer :: x, y
          logical :: b

          b = ASSOCIATED(x % left, y)
        END FUNCTION isLeft

!...........................................................................
        FUNCTION isRight(x, y) result(b)
!..........................................................................
! Check if node y is right child of x
!..........................................................................
          type (treenode), pointer :: x, y
          logical :: b

          b = ASSOCIATED(x % right, y)
        END FUNCTION isRight

!...........................................................................
        FUNCTION isLeaf(x) result(b)
!..........................................................................
! Check if node x is a leaf
!..........................................................................
          type (treenode), pointer :: x
          logical :: b

          b = ASSOCIATED(x, nil)
        END FUNCTION isLeaf
!..........................................................................

        SUBROUTINE tree_collect(this, Nsuper, Nfound, eggs, momfish)
          type (RedBlackTree) :: this
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
          IF (this % nitems .eq. 0) THEN
            RETURN
          ELSE IF (Nsuper >= this % nitems) THEN
            IF (ASSOCIATED(this % root)) CALL copy_traverse(            &
     &               this % root % left, Nsuper, Nfound, eggs, momfish)
            IF (Nfound .ne. this % nitems) print *, "Oh, noes!",        &
     &               Nfound, this % nitems
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
            roe_list % next % egg_tree => this % root
            roe_list % next % egg_sum = this % root % egg_sum
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
              l_tree % right => nil
              l_tree % parent => nil
              r_tree % right => nil
              r_tree % parent => nil

              IF (.not. isLeaf(t_tree % left)) THEN
                l_tree % left => t_tree % left
                t_tree % left => nil
                l_tree % egg_sum = l_tree % left % egg_sum
              ELSE
                l_tree % left => t_tree
                single = .true.
                t_tree % egg_sum = t_tree % eggs
                l_tree % egg_sum = t_tree % eggs
              END IF
              IF (.not. isLeaf(t_tree % right)) THEN
                r_tree % left => t_tree % right
                t_tree % right => nil
                r_tree % egg_sum = r_tree % left % egg_sum
              ELSE
                r_tree % left => t_tree
                single = .true.
                t_tree % egg_sum = t_tree % eggs
                r_tree % egg_sum = t_tree % eggs
              END IF
              l_tree % left % parent => nil
              l_tree % left % red = .false.
              r_tree % left % parent => nil
              r_tree % left % red = .false.
! If we've both subtrees, give the top's eggs to the smaller
! (by egg count) tree, making a new treenode (because that's what
! insert_priv does).
              IF (.not. single) THEN
                IF (l_tree % left % egg_sum <                           &
     &                     r_tree % left % egg_sum) THEN
                  CALL insert_priv(l_tree, t_tree % eggs,               &
     &                 t_tree % dist, t_tree % momfish)
                ELSE
                  CALL insert_priv(r_tree, t_tree % eggs,               &
     &                 t_tree % dist,  t_tree % momfish)
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
        END SUBROUTINE tree_collect

!  Extract the tree with most eggs from the list.
!  The splittable argument if true means that I want the first tree that
!  can be split into sub-trees.
!
!  The objects in egg_batch are the tops of the trees; we are going to
!  return the left child of the top, which is the top node containing
!  actual eggs.
        FUNCTION pop_tree(splittable)
          logical :: splittable
          type(treenode), pointer :: pop_tree
          type(egg_batch), pointer :: cur, prev

          prev => roe_list
          cur => roe_list % next
          DO
!  I'm either returning the first no matter what, or I'm returning the
!  first with more than one node.
            IF (.not. splittable .or. (cur % egg_sum .ne.               &
     &               cur % egg_tree % left % eggs)) THEN
              IF (ASSOCIATED(cur % next)) THEN
                prev % next => cur % next
              END IF
              pop_tree => cur % egg_tree % left
              RETURN
            END IF
            prev => cur
            IF (ASSOCIATED(cur % next)) THEN
              cur => cur % next
            ELSE
              pop_tree => null()
            END IF
          END DO
        END FUNCTION pop_tree

!  Have to keep this sorted by egg count
!  This is going to be easier in C++ with STL.
        SUBROUTINE push_trees(l_tree, r_tree)
          type(treenode), pointer :: l_tree
          type(treenode), pointer :: r_tree
          type(egg_batch), pointer :: left, right
          ALLOCATE(left)
          left % egg_tree => l_tree
          left % egg_sum = l_tree % egg_sum
          CALL push_one_tree(left)
          ALLOCATE(right)
          right % egg_tree => r_tree
          right % egg_sum = r_tree % egg_sum
          CALL push_one_tree(right)
        END SUBROUTINE push_trees

        SUBROUTINE push_one_tree(t_tree)
          type(egg_batch), pointer :: t_tree
          type(egg_batch), pointer :: cur, next
          real(r8) :: cur_eggs, next_eggs, my_eggs

          cur => roe_list
          next => cur % next
          my_eggs = t_tree % egg_sum
          cur_eggs = 1.e35
          IF (ASSOCIATED(next)) THEN
            next_eggs = next % egg_sum
          ELSE
            next_eggs = 0
          END IF
          DO
            IF (my_eggs < cur_eggs .and. my_eggs >= next_eggs) THEN
              t_tree % next => next
              cur % next => t_tree
              RETURN
            END IF
            cur => next
            next => cur % next
            cur_eggs = cur % egg_sum
            IF (ASSOCIATED(next)) THEN
              next_eggs = next % egg_sum
            ELSE
              next_eggs = 0
            END IF
          END DO
        END SUBROUTINE push_one_tree

        SUBROUTINE tree_traverse(this)
          type (RedBlackTree) :: this
          IF (.not. isleaf(this % root)) CALL                           &
     &              traverse_print(this % root % left)
        END SUBROUTINE tree_traverse

        RECURSIVE SUBROUTINE copy_traverse(p, Nsuper, Nfound, eggs,     &
     &                   momfish)
          type(treenode), pointer :: p
          integer, intent(in)  :: Nsuper
          integer, intent(inout) :: Nfound
          real(r8), intent(inout) :: eggs(Nsuper)
          integer, intent(inout) :: momfish(Nsuper)

          IF (.not. isLeaf(p % left)) CALL copy_traverse(p % left,        &
     &                Nsuper, Nfound, eggs, momfish)
          Nfound = Nfound + 1
          eggs(Nfound) = p % eggs
          momfish(Nfound) = p % momfish
          IF (.not. isLeaf(p % right)) CALL copy_traverse(p % right,      &
     &                Nsuper, Nfound, eggs, momfish)
        END SUBROUTINE copy_traverse

        RECURSIVE SUBROUTINE traverse_print(p)
          type(treenode), pointer :: p

          IF (.not. isLeaf(p % left)) CALL traverse_print(p % left)
          CALL print_node(p)
          IF (.not. isLeaf(p % right)) CALL traverse_print(p % right)

        END SUBROUTINE traverse_print

        SUBROUTINE print_node(cur)
          type(treenode), pointer :: cur
          print *, "Node: ", cur % dist, cur % eggs, cur % egg_sum
        END SUBROUTINE print_node

        SUBROUTINE tree_destroy(this)
! Loop through, cleaning up both the treenodes and the listnodes.
          type (RedBlackTree) :: this
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
        END SUBROUTINE tree_destroy

      END MODULE mod_tree
