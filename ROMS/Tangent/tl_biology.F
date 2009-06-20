#include "cppdefs.h"
      MODULE tl_biology_mod
#if defined TANGENT && defined BIOLOGY
!
!svn $Id: tl_biology.F 992 2009-05-29 01:12:37Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine computes the tangent linear source and sink terms for  !
!  the selected biology model.                                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: tl_biology

      CONTAINS

# if defined BIO_FENNEL
!!#  include <tl_fennel.h>
# elif defined NEMURO
!!#  include <tl_nemuro.h>
# elif defined BIO_UMAINE
!!#  include <tl_UMAINE.h>
# elif defined NPZD_FRANKS
#  include <tl_npzd_Franks.h>
# elif defined NPZD_POWELL
#  include <tl_npzd_Powell.h>
# endif

#endif
      END MODULE tl_biology_mod
