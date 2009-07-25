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
!  At the end of each year, we kill off any remaining fish of the
!  the oldest age class to make room for the next age class. We also
!  move the indices for creating new fish.
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

!  Start by killing the old fish for each species.
!  The max age is the same for all species (for now).

      DO isp=1,Nspecies(ng)
        DO i=1,Nfishperyear(ng)
          ifsh = i+if_age_base(isp,Nyearclass(ng)-1)
          FISHES % alive(ifsh) = .FALSE.
          FISHES % bounded(ifsh) = .FALSE.
          DO j=0,NFT
            FISHES(ng) % track(ixgrd,j,ifsh) = -99.
            FISHES(ng) % track(iygrd,j,ifsh) = -99.
            FISHES(ng) % track(izgrd,j,ifsh) = -99.
          END DO
        END DO
      END DO

!  Now we update the indices
      DO isp=1,Nspecies(ng)
         DO i=0,Nyearclass(ng)-1
           if_age_base(isp,i) = if_age_base(isp,i) + Nfishperyear(ng)
           IF (if_age_base(isp,i) .ge. if_sp_base(isp) +                &
     &                  Nfishperyear(ng)*Nyearclass(ng) ) THEN
             if_age_base(isp,i) = if(age_base(isp,i) -                  &
     &                     Nfishperyear(ng)*Nyearclass(ng)
           END IF
         END DO
      END DO

!  Finally, get ready for the new eggs
      DO isp=1,Nspecies(ng)
        FISHES(ng) % num_free(isp) = Nfishperyear(ng)
        FISHES(ng) % num_super(isp) = 0
        FISHES(ng) % next_free(isp) = 1 + if_age_base(isp,0)
      END DO

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE new_year
!
#endif
      END MODULE new_year_mod
