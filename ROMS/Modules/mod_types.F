#include "cppdefs.h"
      MODULE mod_types
!
!svn $Id: mod_ocean.F 975 2009-05-05 22:51:13Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!  Set up data structure for linked lists.                             !
!=======================================================================
!
        implicit none

#if defined NEMURO_SAN
        type fishnode
          type(fishnode), pointer :: next => null()
          integer :: fish
        end type fishnode
# if defined EGGS_BISECTION
        type egg_box
          integer :: i1, i2, j1, j2
          integer :: egg_sum
        end type egg_box
# endif
#endif
      END MODULE mod_types
