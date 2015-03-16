#include "cppdefs.h"
      MODULE mod_types
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
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
# if defined PREDATOR
        type prednode
          type(prednode), pointer :: next => null()
          integer :: pred
        end type prednode
# endif
#endif
      END MODULE mod_types
