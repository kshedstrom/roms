#include "cppdefs.h"
      MODULE new_year_mod
#if defined NEMURO_SAN
!
!svn $Id: step_floats.F 997 2009-06-04 17:42:20Z kate $
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Do what needs to be done at the end of a year to set up for new     !
!  spawning.                                                           !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: new_year

      CONTAINS
!
!***********************************************************************
      SUBROUTINE new_year(ng)
!***********************************************************************
!
      USE mod_param
      USE mod_fish
      USE mod_scalars
      USE mod_stepping
      USE mod_parallel
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_collect
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
      integer :: i, isp, isp2, j, l, Nfound, Navail, Npts, ifsh
      integer :: count(Nspecies(ng)), count_sum
# ifdef DISTRIBUTE
      real(r8) :: Xstr, Xend, Ystr, Yend
      real(r8), dimension(Nfloats(ng)*NFV(ng)*(NFT+1)) :: Fwrk
      real(r8), parameter :: Fspv = 0.0_r8
# endif
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif


# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE new_year
!
#endif
      END MODULE new_year_mod
