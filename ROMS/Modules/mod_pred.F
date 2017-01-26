#include "cppdefs.h"
      MODULE mod_pred
#if defined NEMURO_SAN && defined PREDATOR
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
        USE mod_types
!
        implicit none

        TYPE T_PREDS

          logical, pointer  :: bounded(:)

          integer, pointer  :: Findex(:)
          integer, pointer  :: Ftype(:)

          real(r8), pointer :: Flon(:)
          real(r8), pointer :: Flat(:)
          real(r8), pointer :: Tinfo(:,:)
          real(r8), pointer :: rwalk(:)
          real(r8), pointer :: track(:,:,:)
          real(r8), pointer :: Fz0(:)

! Predator fields
          real(r8), pointer :: bioenergy(:,:)
          integer, pointer  :: species(:)
          logical, pointer  :: alive(:)
          type(prednode), pointer :: prednodes(:)

        END TYPE T_PREDS

        TYPE (T_PREDS), allocatable :: PREDS(:)

      CONTAINS

      SUBROUTINE allocate_pred (ng)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initialize all variables in the module   !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_scalars, only: LdefPRED
      USE mod_biology
      USE mod_fish
      USE mod_floats
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
      IF (ng.eq.1) allocate ( PREDS(Ngrids) )
!
      allocate ( PREDS(ng) % bounded(Npred(ng)) )

      allocate ( PREDS(ng) % Findex(0:Npred(ng)) )

      allocate ( PREDS(ng) % Ftype(Npred(ng)) )

      allocate ( PREDS(ng) % Flon(Npred(ng)) )

      allocate ( PREDS(ng) % Flat(Npred(ng)) )

      allocate ( PREDS(ng) % Tinfo(0:izrhs,Npred(ng)) )

      allocate ( PREDS(ng) % rwalk(Npred(ng)*6) )

      allocate ( PREDS(ng) % track(NFV(ng),0:NFT,Npred(ng)) )

      allocate ( PREDS(ng) % Fz0(Npred(ng)) )

      allocate ( PREDS(ng) % bioenergy(NPredV(ng),Npred(ng)) )

      allocate ( PREDS(ng) % species(Npred(ng)) )

      allocate ( PREDS(ng) % alive(Npred(ng)) )

      allocate ( PREDS(ng) % prednodes(Npred(ng)) )
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      PREDS(ng) % Findex(0) = 0
      DO iflt=1,Npred(ng)
        PREDS(ng) % bounded(iflt) = .FALSE.
        PREDS(ng) % Findex(iflt) = 0
        PREDS(ng) % Ftype(iflt) = 0
        PREDS(ng) % Flon(iflt) = IniVal
        PREDS(ng) % Flat(iflt) = IniVal
        PREDS(ng) % Fz0(iflt) = 0
        PREDS(ng) % rwalk(iflt) = IniVal
        PREDS(ng) % rwalk(iflt+Npred(ng)) = IniVal
        PREDS(ng) % rwalk(iflt+2*Npred(ng)) = IniVal
        PREDS(ng) % rwalk(iflt+3*Npred(ng)) = IniVal
        PREDS(ng) % rwalk(iflt+4*Npred(ng)) = IniVal
        PREDS(ng) % rwalk(iflt+5*Npred(ng)) = IniVal
        DO i=0,izrhs
          PREDS(ng) % Tinfo(i,iflt) = IniVal
        END DO
        DO i=1,NFV(ng)
          PREDS(ng) % track(i,0,iflt) = IniVal
          PREDS(ng) % track(i,1,iflt) = IniVal
          PREDS(ng) % track(i,2,iflt) = IniVal
          PREDS(ng) % track(i,3,iflt) = IniVal
          PREDS(ng) % track(i,4,iflt) = IniVal
        END DO
        DO i=1,NPredV(ng)
          PREDS(ng) % bioenergy(i,iflt) = IniVal
        END DO
        PREDS(ng) % species(iflt) = if_none
        PREDS(ng) % alive(iflt) = .FALSE.
        PREDS(ng) % prednodes(iflt) % pred = iflt
      END DO

      LdefPRED(ng)=.TRUE.

      RETURN
      END SUBROUTINE allocate_pred
#endif
      END MODULE mod_pred
