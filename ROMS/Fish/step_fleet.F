#include "cppdefs.h"
      MODULE step_fleet_mod
#if defined NEMURO_SAN && defined FISHING_FLEET
!
!svn $Id$
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Compute fishing from fleet                                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: step_fleet

      CONTAINS
!
!***********************************************************************
      SUBROUTINE step_fleet(ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_grid
      USE mod_fish
      USE mod_fleet
      USE mod_fishing
      USE mod_biology
      USE mod_scalars
      USE mod_stepping
      USE mod_parallel
      USE nrutil
      USE ran_state, ONLY: ran_seed
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: ib, ip, ifid
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif
!
      CALL gather_array(ng, iNLM, LBi, UBi, LBj, UBj,                   &
     &                 GRID(ng)%xr, GRID(ng)%yr, OCEAN(ng)%fish_count,  &
     &                 OCEAN(ng)%distport, OCEAN(ng)%nearport,          &
     &                 OCEAN(ng)%cpue, OCEAN(ng)%avgcpue,               &
                       OCEAN(ng)%percpue, OCEAN(ng)%new_climcpue,       &
                       OCEAN(ng)%climcpue)
!
      IF (Master) THEN
!
!----  Intialization of port distances ----
        IF (BOATS(ng)%initports) THEN
          CALL initports(ng)
          BOATS(ng)%initports=.FALSE.
        END IF
!
!----  Intialization of CPUE based on fish-by-cell ----
       IF (BOATS(ng)%initcpue) THEN
         CALL initcpue(ng)
         BOATS(ng)%initcpue=.FALSE.
       END IF
!
!----  Reset fish mortality ----
        DO ifid=1,Nfish(ng)
          FISHES(ng)%fmortF(ifid)=0.0_r8
        END DO

!----  Compute daily catch for each boat ----
        DO ib=1,Nboats(ng)
! Reset time counter to 24h
          BOATS(ng)%tAvail=24.0_r8
          BOATS(ng)%catch(ib)=0.0_r8
          BOATS(ng)%boat(ibnloc,ib)=0
          BOATS(ng)%boat(ibifsh1,ib)=0
          BOATS(ng)%boat(ibjfsh1,ib)=0
          BOATS(ng)%boat(ibifsh2,ib)=0
          BOATS(ng)%boat(ibjfsh2,ib)=0
          CALL bestlocation(ng,ib)
! At least on cell has a positive net revenue based on expected CPUE
          IF (BOATS(ng)%blocflag(ib).eq.1) THEN
            BOATS(ng)%boat(ibnloc,ib)=BOATS(ng)%boat(ibnloc,ib)+1
            BOATS(ng)%boat(ibifsh1,ib)=BOATS(ng)%boat(ibiloc,ib)
            BOATS(ng)%boat(ibjfsh1,ib)=BOATS(ng)%boat(ibjloc,ib)
            CALL catch(ng,ib)
! No fish in cell or sardine encountered -> Try a 2nd location
            IF ((BOATS(ng)%bfishflag(ib).eq.0).or.                      &
     &          (BOATS(ng)%bfishflag(ib).eq.1)) THEN
              CALL bestlocation(ng,ib)
              BOATS(ng)%boat(ibnloc,ib)=BOATS(ng)%boat(ibnloc,ib)+1
              BOATS(ng)%boat(ibifsh2,ib)=BOATS(ng)%boat(ibiloc,ib)
              BOATS(ng)%boat(ibjfsh2,ib)=BOATS(ng)%boat(ibjloc,ib)
              CALL catch(ng,ib)
            END IF
! No time or spcae left -> Done fishing
            IF (BOATS(ng)%bfishflag(ib).eq.2) THEN
              CONTINUE
            END IF
! No more sardine in cell, but time and space left -> Try a 2nd location
            IF (BOATS(ng)%bfishflag(ib).eq.3) THEN
              CALL bestlocation(ng,ib)
! At least on cell has a positive net revenue based on expected CPUE
              IF (BOATS(ng)%blocflag(ib).eq.1) THEN
                BOATS(ng)%boat(ibnloc,ib)=BOATS(ng)%boat(ibnloc,ib)+1
                BOATS(ng)%boat(ibifsh2,ib)=BOATS(ng)%boat(ibiloc,ib)
                BOATS(ng)%boat(ibjfsh2,ib)=BOATS(ng)%boat(ibjloc,ib)
                CALL catch(ng,ib)
              END IF
            END IF
          END IF
! Go to nearest port (nearest port id already set in "bestlocation")
          ip=BOATS(ng)%boat(ibport,ib)
          BOATS(ng)%boat(ibiloc,ib)=iPort(ip,ng)
          BOATS(ng)%boat(ibjloc,ib)=jPort(ip,ng)
!        print*, 'CATCH', ib, BOATS(ng)%catch(ib), BOATS(ng)%tAvail(ib)
        END DO
!
!----   Aggregated catch values  ----
!  Note: TotCatch, TotLoc, PortCatch, and AvgTotCatch can be done as postprocessing
        BOATS(ng)%TotCatch=0.0_r8
        BOATS(ng)%TotLoc=0.0_r8
        DO ip=1,Nports(ng)
          BOATS(ng)%PortCatch(ip)=0.0_r8
        END DO
        DO ib=1,Nboats(ng)
          BOATS(ng)%cumcatch(ib)=BOATS(ng)%cumcatch(ib)+                &
     &                               BOATS(ng)%catch(ib)
          BOATS(ng)%TotCatch=BOATS(ng)%TotCatch+BOATS(ng)%catch(ib)
          BOATS(ng)%TotLoc=BOATS(ng)%TotLoc+BOATS(ng)%boat(ibnloc,ib)
          ip=BOATS(ng)%boat(ibport,ib)
          BOATS(ng)%PortCatch(ip)=BOATS(ng)%PortCatch(ip)+              &
     &                                 BOATS(ng)%catch(ib)
        END DO
        BOATS(ng)%AnnCatch=BOATS(ng)%AnnCatch+BOATS(ng)%TotCatch
        IF (BOATS(ng)%TotLoc.gt.0) THEN
          BOATS(ng)%AvgTotCatch=BOATS(ng)%TotCatch/                     &
     &                          REAL(BOATS(ng)%TotLoc,r8)
        ELSE
          BOATS(ng)%AvgTotCatch=0.0_r8
        END IF
!
        IF(BOATS(ng)%TotFishBio.gt.0.001_r8) THEN
          BOATS(ng)%AnnFishMort=-LOG(1.0_r8-                            &
     &                   BOATS(ng)%AnnCatch/BOATS(ng)%TotFishBio)
        ELSE
          BOATS(ng)%AnnFishMort=0.0_r8
        END IF
!        print*, 'SUMCATCH', BOATS(ng)%TotCatch, BOATS(ng)%TotLoc
!
!----  Update CPUE to reflect today's catch  ----
! (yearly update of climatological CPUE done in new_year_fleet)
        CALL updatecpue(ng)
!
      END IF
!
      CALL scatter_array(ng, iNLM, LBi, UBi, LBj, UBj,                  &
     &                   OCEAN(ng)%distport, OCEAN(ng)%nearport,        &
     &                   OCEAN(ng)%cpue, OCEAN(ng)%avgcpue,             &
     &                   OCEAN(ng)%percpue, OCEAN(ng)%new_climcpue,     &
     &                   OCEAN(ng)%climcpue)
!
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE step_fleet
!
#endif
      END MODULE step_fleet_mod
