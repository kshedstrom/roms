#include "cppdefs.h"
      MODULE biology_mod
#if defined NONLINEAR && defined BIOLOGY
!
!svn $Id: biology.F 984 2009-05-24 01:43:19Z kate $
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

# if defined BIO_FENNEL
#  include <fennel.h>
# elif defined NEMURO
#  include <nemuro.h>
# elif defined NPZD_FRANKS
#  include <npzd_Franks.h>
# elif defined NPZD_POWELL
#  include <npzd_Powell.h>
# elif defined ECOSIM
#  include <ecosim.h>
# elif defined BIO_UMaine
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
