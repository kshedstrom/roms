#include "cppdefs.h"
#ifdef STOPERTURB
      MODULE mod_stoperturb
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!                   Stochastic perturbations                           !
!                                                                      !
!  Forcing fields.                                                     !
!                                                                      !
!  Uwind        Surface wind in the XI-direction (m/s) at              !
!                 horizontal RHO-points.                               !
!  UwindG       Latest two-time snapshots of input "Uwind" grided      !
!                 data used for interpolation.                         !
!  Vwind        Surface wind in the ETA-direction (m/s) at             !
!                 horizontal RHO-points.                               !
!  VwindG       Latest two-time snapshots of input "Vwind" grided      !
!                 data used for interpolation.                         !
!  Tair         Surface air temperature (Celsius)                      !
!  TairG        Latest two-time snapshots of input "Tair" grided       !
!                 data used for interpolation.                         !
!  srflx        Kinematic surface shortwave solar radiation flux       !
!                 (Celsius m/s) at horizontal RHO-points               !
!  srflxG       Latest two-time snapshots of input "srflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!=======================================================================
!

        USE mod_kinds
        USE mod_scalars
        implicit none

        TYPE T_PERTURB

#ifdef FORCE_PERTURB
          real(r8), pointer :: UwindG(:,:,:)
          real(r8), pointer :: VwindG(:,:,:)
          real(r8), pointer :: TairG(:,:,:)
          real(r8), pointer :: srflxG(:,:,:)
#endif

#ifdef EOS_PERTURB
          real(r8), pointer :: deltax(:,:)
          real(r8), pointer :: deltay(:,:)
          real(r8), pointer :: deltaz(:,:)
#endif

        END TYPE T_PERTURB

        TYPE (T_PERTURB),     allocatable :: PERTURB(:)

        integer, allocatable :: iDforcperturb(:)
        integer              :: Nforcperturb=4

      CONTAINS


      SUBROUTINE allocate_perturb (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) THEN
         allocate ( PERTURB(Ngrids) )
      ENDIF
!
!  Nonlinear model state
!
#ifdef FORCE_PERTURB
      allocate ( PERTURB(ng) % UwindG(LBi:UBi,LBj:UBj,2) )
      allocate ( PERTURB(ng) % VwindG(LBi:UBi,LBj:UBj,2) )
      allocate ( PERTURB(ng) % TairG(LBi:UBi,LBj:UBj,2) )
      allocate ( PERTURB(ng) % srflxG(LBi:UBi,LBj:UBj,2) )
      allocate ( iDforcperturb(Nforcperturb) )
#endif
#ifdef EOS_PERTURB
      allocate ( PERTURB(ng) % deltax(LBi:UBi,LBj:UBj) )
      allocate ( PERTURB(ng) % deltay(LBi:UBi,LBj:UBj) )
      allocate ( PERTURB(ng) % deltaz(LBi:UBi,LBj:UBj) )
#endif
      RETURN
      END SUBROUTINE allocate_perturb


      SUBROUTINE initialize_perturb (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE pseudo_rand2D, ONLY : random_seed_fixed
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k
      real(r8), parameter :: IniVal = 0.0_r8

#include "set_bounds.h"


!
!  Initialize seed for random number generation
!
      CALL random_seed_fixed(seed_user(ng))

!
!  Set array initialization range.
!
#ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#else
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
#ifdef FORCE_PERTURB
            PERTURB(ng) % UwindG(i,j,1) = IniVal
            PERTURB(ng) % UwindG(i,j,2) = IniVal
            PERTURB(ng) % VwindG(i,j,1) = IniVal
            PERTURB(ng) % VwindG(i,j,2) = IniVal
            PERTURB(ng) % TairG(i,j,1) = IniVal
            PERTURB(ng) % TairG(i,j,2) = IniVal
            PERTURB(ng) % srflxG(i,j,1) = IniVal
            PERTURB(ng) % srflxG(i,j,2) = IniVal
#endif
#ifdef EOS_PERTURB
            PERTURB(ng) % deltax(i,j) = IniVal
            PERTURB(ng) % deltay(i,j) = IniVal
            PERTURB(ng) % deltaz(i,j) = IniVal
#endif
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE initialize_perturb

      END MODULE mod_stoperturb
#endif

