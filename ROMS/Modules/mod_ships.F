#include "cppdefs.h"
      MODULE mod_ships
#ifdef NEMURO_SAN
!
!svn $Id: mod_fish.F 975 2009-05-05 22:51:13Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
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
!  track      Multivariate float trajectory data at several time       !
!               time levels.                                           !
!                                                                      !
!=======================================================================
!
        USE mod_param
        USE mod_types
!
        implicit none

        TYPE T_SHIPS

          logical, pointer  :: bounded(:)

          integer, pointer  :: Findex(:)
          integer, pointer  :: Ftype(:)

          real(r8), pointer :: Flon(:)
          real(r8), pointer :: Flat(:)
          real(r8), pointer :: Tinfo(:,:)
          real(r8), pointer :: track(:,:,:)
          real(r8), pointer :: Fz0(:)

! Ship fields
!         real(r8), pointer :: bioenergy(:,:)
          integer, pointer  :: species(:)
	  type(fishnode), pointer :: shipnodes(:)

        END TYPE T_SHIPS

        TYPE (T_SHIPS), allocatable :: SHIPS(:)

      CONTAINS

      SUBROUTINE allocate_ships (ng)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initialize all variables in the module   !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
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
      IF (ng.eq.1) allocate ( SHIPS(Ngrids) )
!
      allocate ( SHIPS(ng) % bounded(Nships(ng)) )

      allocate ( SHIPS(ng) % Findex(0:Nships(ng)) )

      allocate ( SHIPS(ng) % Ftype(Nships(ng)) )

      allocate ( SHIPS(ng) % Flon(Nships(ng)) )

      allocate ( SHIPS(ng) % Flat(Nships(ng)) )

      allocate ( SHIPS(ng) % Tinfo(0:izrhs,Nships(ng)) )

      allocate ( SHIPS(ng) % track(NFV(ng),0:NFT,Nships(ng)) )

      allocate ( SHIPS(ng) % Fz0(Nships(ng)) )

!     allocate ( SHIPS(ng) % feedback(NBT,Nships(ng)) )

!     allocate ( SHIPS(ng) % bioenergy(NFishV(ng),Nships(ng)) )

!     allocate ( SHIPS(ng) % species(Nships(ng)) )

!     allocate ( SHIPS(ng) % lifestage(Nships(ng)) )

!     allocate ( SHIPS(ng) % alive(Nships(ng)) )

      allocate ( SHIPS(ng) % shipnodes(Nships(ng)) )
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      SHIPS(ng) % Findex(0) = 0
      DO iflt=1,Nships(ng)
        SHIPS(ng) % bounded(iflt) = .FALSE.
        SHIPS(ng) % Findex(iflt) = 0
        SHIPS(ng) % Ftype(iflt) = 0
        SHIPS(ng) % Flon(iflt) = IniVal
        SHIPS(ng) % Flat(iflt) = IniVal
        SHIPS(ng) % Fz0(iflt) = 0
!       DO i=1,NBT
!         SHIPS(ng) % feedback(i,iflt) = IniVal
!       END DO
        DO i=0,izrhs
          SHIPS(ng) % Tinfo(i,iflt) = IniVal
        END DO
        DO i=1,NFV(ng)
          SHIPS(ng) % track(i,0,iflt) = IniVal
          SHIPS(ng) % track(i,1,iflt) = IniVal
          SHIPS(ng) % track(i,2,iflt) = IniVal
          SHIPS(ng) % track(i,3,iflt) = IniVal
          SHIPS(ng) % track(i,4,iflt) = IniVal
        END DO
!       DO i=1,NFishV(ng)
!         SHIPS(ng) % bioenergy(i,iflt) = IniVal
!       END DO
!       SHIPS(ng) % species(iflt) = if_none
!       SHIPS(ng) % lifestage(iflt) = if_none
!       SHIPS(ng) % alive(iflt) = .FALSE.
        SHIPS(ng) % shipnodes(iflt) % fish = iflt
      END DO

      RETURN
      END SUBROUTINE allocate_ships
#endif
      END MODULE mod_ships
