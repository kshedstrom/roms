#include "cppdefs.h"
      MODULE mod_fleet
#if defined NEMURO_SAN && defined FISHING_FLEET
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Findex     Indices of spherical coordinates entries in initial      !
!               location arrays, if any.                               !
!  Flon       Initial longitude locations, if any.                     !
!  Flat       Initial latitude locations, if any.                      !
!  Ftype      Float trajectory type:                                   !
!               Ftype(:) = 1,  neutral density 3D Lagrangian           !
!               Ftype(:) = 2,  isobaric (constant depth) float.        !
!  Tinfo      Float trajectory initial information.                    !
!  bounded    Float bounded status switch.                             !
!  rwalk      Normally distributed random deviates used in 3-D         !
!               random walk.                                           !
!  track      Multivariate float trajectory data at several time       !
!               time levels.                                           !
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none

        TYPE T_BOATS

! Boat fields
          integer, pointer :: boat(:,:)
          integer, pointer :: blocflag(:)
          integer, pointer :: bfishflag(:)
          integer, pointer :: TotLoc
          logical, pointer :: initcpue
          logical, pointer :: initports
          real(r8), pointer :: tAvail(:)
          real(r8), pointer :: t2LocOpt(:)
          real(r8), pointer :: t2PortOpt(:)
          real(r8), pointer :: catch(:)
          real(r8), pointer :: cumcatch(:)
          real(r8), pointer :: TotCatch
          real(r8), pointer :: AvgTotCatch
          real(r8), pointer :: AnnCatch
          real(r8), pointer :: TotFishBio
          real(r8), pointer :: AnnFishMort
          real(r8), pointer :: PortCatch(:)
          real(r8), pointer :: rwalk(:)

        END TYPE T_BOATS

        TYPE (T_BOATS), allocatable :: BOATS(:)

      CONTAINS

      SUBROUTINE allocate_fleet (ng)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initialize all variables in the module   !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
      USE mod_biology
      USE mod_fish
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: i, iflt

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( BOATS(Ngrids) )
!
      allocate ( BOATS(ng) % boat(NBoatV(ng),Nboats(ng)) )

      allocate ( BOATS(ng) % blocflag(Nboats(ng)) )

      allocate ( BOATS(ng) % bfishflag(Nboats(ng)) )

      allocate ( BOATS(ng) % TotLoc )

      allocate ( BOATS(ng) % initcpue )

      allocate ( BOATS(ng) % initports )

      allocate ( BOATS(ng) % tAvail(Nboats(ng)) )

      allocate ( BOATS(ng) % t2LocOpt(Nboats(ng)) )

      allocate ( BOATS(ng) % t2PortOpt(Nboats(ng)) )

      allocate ( BOATS(ng) % catch(Nboats(ng)) )

      allocate ( BOATS(ng) % cumcatch(Nboats(ng)) )

      allocate ( BOATS(ng) % TotCatch )

      allocate ( BOATS(ng) % AvgTotCatch )

      allocate ( BOATS(ng) % AnnCatch )

      allocate ( BOATS(ng) % TotFishBio )

      allocate ( BOATS(ng) % AnnFishMort )

      allocate ( BOATS(ng) % PortCatch(Nports(ng)) )

      allocate ( BOATS(ng) % rwalk(Nboats(ng)) )
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO iflt=1,Nboats(ng)
        DO i=1,NBoatV(ng)
          BOATS(ng) % boat(i,iflt) = IniVal
        END DO
        BOATS(ng) % blocflag(iflt) = 0
        BOATS(ng) % bfishflag(iflt) = 0
        BOATS(ng) % tAvail(iflt) = IniVal
        BOATS(ng) % t2LocOpt(iflt) = IniVal
        BOATS(ng) % t2PortOpt(iflt) = IniVal
        BOATS(ng) % catch(iflt) = IniVal
        BOATS(ng) % cumcatch(iflt) = IniVal
        BOATS(ng) % rwalk(iflt) = IniVal
      END DO
      BOATS(ng) % initcpue = .FALSE.
      BOATS(ng) % initports = .FALSE.
      BOATS(ng) % TotLoc = 0
      BOATS(ng) % TotCatch = IniVal
      BOATS(ng) % AvgTotCatch = IniVal
      BOATS(ng) % AnnCatch = IniVal
      BOATS(ng) % TotFishBio = IniVal
      BOATS(ng) % AnnFishMort = IniVal
      DO iflt=1,Nports(ng)
        BOATS(ng) % PortCatch(iflt) = IniVal
      END DO

!RD: add this
      LdefFLEET(ng)=.TRUE.


      RETURN
      END SUBROUTINE allocate_fleet
#endif
      END MODULE mod_fleet
