#include "cppdefs.h"
      MODULE biology_mod
#if defined NONLINEAR && defined BIOLOGY
!
!svn $Id: biology.F 1076 2009-09-25 23:18:43Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the source and sink terms for selected        !
!   biology model.                                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: biology

      CONTAINS
# if defined BEST_NPZ
#  include <bestnpz.h>
# elif defined BIO_FENNEL
#  include <fennel.h>
# elif defined ECOSIM
#  include <ecosim.h>
# elif defined NEMURO
#  include <nemuro.h>
# elif defined NPZD_FRANKS
#  include <npzd_Franks.h>
# elif defined NPZD_IRON
#  include <npzd_iron.h>
# elif defined NPZD_POWELL
#  include <npzd_Powell.h>
# elif defined BIO_UMAINE
#  include <bio_UMAINE.h>
# elif defined BIO_GOANPZ
#  include <goanpz.h>
# else
      SUBROUTINE biology(ng, tile)
      integer ng, tile
      END subroutine biology
# endif

#endif
      END MODULE biology_mod
