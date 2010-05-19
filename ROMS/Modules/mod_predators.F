#include "cppdefs.h"
      MODULE mod_predators
#ifdef NEMURO_SAN
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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

        TYPE T_PREDATORS

          logical, pointer  :: bounded(:)

          integer, pointer  :: Findex(:)
          integer, pointer  :: Ftype(:)

          real(r8), pointer :: Flon(:)
          real(r8), pointer :: Flat(:)
          real(r8), pointer :: Tinfo(:,:)
          real(r8), pointer :: track(:,:,:)
          real(r8), pointer :: Fz0(:)

! Predator fields
          real(r8), pointer :: bioenergy(:,:)
          real(r8), pointer :: feedback(:,:)
          integer, pointer  :: species(:)
          logical, pointer  :: alive(:)
          logical, pointer  :: eatfish(:)
          type(fishnode), pointer :: fishnodes(:)

        END TYPE T_PREDATORS

        TYPE (T_PREDATORS), allocatable :: PREDATORS(:)

      CONTAINS

      SUBROUTINE allocate_predators (ng)
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
      IF (ng.eq.1) allocate ( PREDATORS(Ngrids) )
!
      allocate ( PREDATORS(ng) % bounded(Npred(ng)) )

      allocate ( PREDATORS(ng) % Findex(0:Npred(ng)) )

      allocate ( PREDATORS(ng) % Ftype(Npred(ng)) )

      allocate ( PREDATORS(ng) % Flon(Npred(ng)) )

      allocate ( PREDATORS(ng) % Flat(Npred(ng)) )

      allocate ( PREDATORS(ng) % Tinfo(0:izrhs,Npred(ng)) )

      allocate ( PREDATORS(ng) % track(NFV(ng),0:NFT,Npred(ng)) )

      allocate ( PREDATORS(ng) % Fz0(Npred(ng)) )

      allocate ( PREDATORS(ng) % feedback(NT(ng),Npred(ng)) )

      allocate ( PREDATORS(ng) % bioenergy(NFishV(ng),Npred(ng)) )

      allocate ( PREDATORS(ng) % species(Npred(ng)) )

      allocate ( PREDATORS(ng) % alive(Npred(ng)) )

      allocate ( PREDATORS(ng) % eatfish(Npred(ng)) )

      allocate ( PREDATORS(ng) % fishnodes(Npred(ng)) )
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      PREDATORS(ng) % Findex(0) = 0
      DO iflt=1,Npred(ng)
        PREDATORS(ng) % bounded(iflt) = .FALSE.
        PREDATORS(ng) % Findex(iflt) = 0
        PREDATORS(ng) % Ftype(iflt) = 0
        PREDATORS(ng) % Flon(iflt) = IniVal
        PREDATORS(ng) % Flat(iflt) = IniVal
        PREDATORS(ng) % Fz0(iflt) = 0
        DO i=1,NT(ng)
          PREDATORS(ng) % feedback(i,iflt) = IniVal
        END DO
        DO i=0,izrhs
          PREDATORS(ng) % Tinfo(i,iflt) = IniVal
        END DO
        DO i=1,NFV(ng)
          PREDATORS(ng) % track(i,0,iflt) = IniVal
          PREDATORS(ng) % track(i,1,iflt) = IniVal
          PREDATORS(ng) % track(i,2,iflt) = IniVal
          PREDATORS(ng) % track(i,3,iflt) = IniVal
          PREDATORS(ng) % track(i,4,iflt) = IniVal
        END DO
        DO i=1,NFishV(ng)
          PREDATORS(ng) % bioenergy(i,iflt) = IniVal
        END DO
        PREDATORS(ng) % species(iflt) = if_none
        PREDATORS(ng) % alive(iflt) = .FALSE.
        PREDATORS(ng) % eatfish(iflt) = .FALSE.
        PREDATORS(ng) % fishnodes(iflt) % fish = iflt
      END DO

      RETURN
      END SUBROUTINE allocate_predators
#endif
      END MODULE mod_predators
