#include "cppdefs.h"
      MODULE new_year_mod
#if defined NEMURO_SAN
!
!svn $Id$
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
      real(r8), parameter :: Fspv = 0.0_r8
# endif
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif

!  Update the indices
      DO isp=1,Nspecies(ng)
        FISHES(ng)%age_base(isp)=FISHES(ng)%age_base(isp)+              &
     &                                        Nfishperyear(ng)
        IF (FISHES(ng)%age_base(isp).ge.FISHES(ng)%species_base(isp)+   &
     &                            Nfishperyear(ng)*Nyearclass(ng)) THEN
          FISHES(ng)%age_base(isp)=FISHES(ng)%species_base(isp)
        END IF
      END DO

!  Killing the old fish for each species.
!  The max age is the same for all species (for now).

      DO isp=1,Nspecies(ng)
        DO i=1,Nfishperyear(ng)
          ifsh = i+FISHES(ng) % age_base(isp)
          FISHES(ng) % alive(ifsh) = .FALSE.
          FISHES(ng) % bounded(ifsh) = .FALSE.
          FISHES(ng) % bioenergy(ifwwt,ifsh) = 0.0_r8
          FISHES(ng) % bioenergy(ifworth,ifsh) = 0.0_r8
          FISHES(ng) % bioenergy(ifage,ifsh) = 0.0_r8
          FISHES(ng) % bioenergy(ifbday,ifsh) = 0.0_r8
          FISHES(ng) % lifestage(ifsh) = 0
          DO j=0,NFT
            FISHES(ng) % track(ixgrd,j,ifsh) = 0.
            FISHES(ng) % track(iygrd,j,ifsh) = 0.
            FISHES(ng) % track(izgrd,j,ifsh) = 0.
          END DO
        END DO
      END DO

!  Get ready for the new eggs
      DO isp=1,Nspecies(ng)
        FISHES(ng)%num_free(isp)=Nfishperyear(ng)
        FISHES(ng)%num_super(isp)=0
        FISHES(ng)%next_free(isp)=1+FISHES(ng)%age_base(isp)
      END DO
      IF (Master) THEN
        print *, 'New Year for Fish'
	print *, 'Free fish ', FISHES(ng) % num_free
	print *, 'Fish ages ', FISHES(ng) % age_base
      END IF

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE new_year
!
#endif
      END MODULE new_year_mod
