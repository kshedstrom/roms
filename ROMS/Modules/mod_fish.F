#include "cppdefs.h"
      MODULE mod_fish
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

        TYPE T_FISHES

          logical, pointer  :: bounded(:)

          integer, pointer  :: Findex(:)
          integer, pointer  :: Ftype(:)

          real(r8), pointer :: Flon(:)
          real(r8), pointer :: Flat(:)
          real(r8), pointer :: Tinfo(:,:)
          real(r8), pointer :: track(:,:,:)
          real(r8), pointer :: Fz0(:)

! Fish fields
          real(r8), pointer :: bioenergy(:,:)
          real(r8), pointer :: feedback(:,:)
          integer, pointer  :: species(:)
          integer, pointer  :: lifestage(:)
          integer, pointer  :: swimtype(:)
          logical, pointer  :: alive(:)
          logical, pointer  :: eatfish(:)
          type(fishnode), pointer :: fishnodes(:)
          real(r8), pointer :: sp_start(:)
          real(r8), pointer :: sp_end(:)
          integer, pointer  :: num_free(:)
          integer, pointer  :: next_free(:)
          integer, pointer  :: num_super(:)

        END TYPE T_FISHES

        TYPE (T_FISHES), allocatable :: FISHES(:)

      CONTAINS

      SUBROUTINE allocate_fish (ng)
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
      IF (ng.eq.1) allocate ( FISHES(Ngrids) )
!
      allocate ( FISHES(ng) % bounded(Nfish(ng)) )

      allocate ( FISHES(ng) % Findex(0:Nfish(ng)) )

      allocate ( FISHES(ng) % Ftype(Nfish(ng)) )

      allocate ( FISHES(ng) % Flon(Nfish(ng)) )

      allocate ( FISHES(ng) % Flat(Nfish(ng)) )

      allocate ( FISHES(ng) % Tinfo(0:izrhs,Nfish(ng)) )

      allocate ( FISHES(ng) % track(NFV(ng),0:NFT,Nfish(ng)) )

      allocate ( FISHES(ng) % Fz0(Nfish(ng)) )

      allocate ( FISHES(ng) % feedback(NT(ng),Nfish(ng)) )

      allocate ( FISHES(ng) % bioenergy(NFishV(ng),Nfish(ng)) )

      allocate ( FISHES(ng) % species(Nfish(ng)) )

      allocate ( FISHES(ng) % lifestage(Nfish(ng)) )

      allocate ( FISHES(ng) % swimtype(Nfish(ng)) )

      allocate ( FISHES(ng) % alive(Nfish(ng)) )

      allocate ( FISHES(ng) % eatfish(Nfish(ng)) )

      allocate ( FISHES(ng) % fishnodes(Nfish(ng)) )

      allocate ( FISHES(ng) % sp_start(Nspecies(ng)) )

      allocate ( FISHES(ng) % sp_end(Nspecies(ng)) )

      allocate ( FISHES(ng) % num_free(Nspecies(ng)) )

      allocate ( FISHES(ng) % num_super(Nspecies(ng)) )

      allocate ( FISHES(ng) % next_free(Nspecies(ng)) )
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      FISHES(ng) % Findex(0) = 0
      DO iflt=1,Nfish(ng)
        FISHES(ng) % bounded(iflt) = .FALSE.
        FISHES(ng) % Findex(iflt) = 0
        FISHES(ng) % Ftype(iflt) = 0
        FISHES(ng) % Flon(iflt) = IniVal
        FISHES(ng) % Flat(iflt) = IniVal
        FISHES(ng) % Fz0(iflt) = 0
        DO i=1,NT(ng)
          FISHES(ng) % feedback(i,iflt) = IniVal
        END DO
        DO i=0,izrhs
          FISHES(ng) % Tinfo(i,iflt) = IniVal
        END DO
        DO i=1,NFV(ng)
          FISHES(ng) % track(i,0,iflt) = IniVal
          FISHES(ng) % track(i,1,iflt) = IniVal
          FISHES(ng) % track(i,2,iflt) = IniVal
          FISHES(ng) % track(i,3,iflt) = IniVal
          FISHES(ng) % track(i,4,iflt) = IniVal
        END DO
        DO i=1,NFishV(ng)
          FISHES(ng) % bioenergy(i,iflt) = IniVal
        END DO
        FISHES(ng) % species(iflt) = if_none
        FISHES(ng) % lifestage(iflt) = if_none
        FISHES(ng) % swimtype(iflt) = 0
        FISHES(ng) % alive(iflt) = .FALSE.
        FISHES(ng) % eatfish(iflt) = .FALSE.
        FISHES(ng) % fishnodes(iflt) % fish = iflt
      END DO

      RETURN
      END SUBROUTINE allocate_fish
#endif
      END MODULE mod_fish
