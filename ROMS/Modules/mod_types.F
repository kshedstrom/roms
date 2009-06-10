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
          type(fishnode), pointer :: next
          integer :: fish
        end type fishnode

        type treenode
          type(treenode), pointer :: left, right
          logical :: red
          integer :: eggs
        end type treenode
#endif
      END MODULE mod_types
