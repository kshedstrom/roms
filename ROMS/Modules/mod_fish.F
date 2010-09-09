#include "cppdefs.h"
      MODULE mod_fish
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
          logical, pointer  :: fishgirl(:)
          type(fishnode), pointer :: fishnodes(:)
          real(r8), pointer :: egg_dur(:)
          real(r8), pointer :: egg_num(:)
          real(r8), pointer :: ysac_dur(:)
          real(r8), pointer :: ysac_num(:)
          real(r8), pointer :: larv_dur(:)
          real(r8), pointer :: larv_num(:)
          real(r8), pointer :: juv_dur(:)
          real(r8), pointer :: juv_num(:)
          real(r8), pointer :: suba_num(:)
          real(r8), pointer :: sp_start(:)
          real(r8), pointer :: sp_end(:)
          real(r8), pointer :: fmortN(:)
          real(r8), pointer :: fmortS(:)
          real(r8), pointer :: fmortF(:)
          real(r8), pointer :: fmortP(:)
          integer, pointer  :: num_free(:)
          integer, pointer  :: next_free(:)
          integer, pointer  :: num_super(:)
          integer, pointer  :: cellid(:)

! "Pointers" into big Nfish dimension:
          integer, pointer :: age_base(:)
          integer, pointer :: species_base(:)
        END TYPE T_FISHES

        TYPE (T_FISHES), allocatable :: FISHES(:)
!
!  Set biological tracer identification indices.
!
      integer, parameter :: ifwwt = 1          ! fish wet weight
      integer, parameter :: ifworth = 2        ! fish worth
      integer, parameter :: ifage = 3          ! age as stage
      integer, parameter :: ifbday = 4         ! birthday
      integer, parameter :: ifspwnloc = 5      ! spawn_dist of mother
      integer, parameter :: ifbatch = 6        ! number of batches spawned
      integer, parameter :: iftspwn = 7        ! time of last batch spawned
      integer, parameter :: ifeggs = 8         ! number of eggs spawned
      integer, parameter :: ifpval = 9         ! fish p-value
      integer, parameter :: ifcsmPS = 10       ! PS consumption
      integer, parameter :: ifcsmPL = 11       ! PL consumption
      integer, parameter :: ifcsmZS = 12       ! ZS consumption
      integer, parameter :: ifcsmZL = 13       ! ZL consumption
      integer, parameter :: ifcsmZP = 14       ! ZP consumption
      integer, parameter :: ifresp = 15        ! Respiration
! Lifestage and gender
      integer, parameter :: if_egg = 1
      integer, parameter :: if_yolksac = 2
      integer, parameter :: if_larva = 3
      integer, parameter :: if_juvenile = 4
      integer, parameter :: if_subadult = 5
      integer, parameter :: if_adult = 6
! Species IDs
      integer, parameter :: if_none = 0
      integer, parameter :: if_anchovy = 1
      integer, parameter :: if_sardine = 2
      integer, parameter :: if_herring = 3
      integer, parameter :: if_polluck = 4
      integer, parameter :: if_pinksalmon = 5
      integer, parameter :: max_yearclass = 10
# ifdef PREDATOR
      integer, parameter :: ipwwt = 1         ! pred wet weight
      integer, parameter :: ipworth = 2       ! pred worth
      integer, parameter :: ippval = 3        ! pred p-value
      integer, parameter :: ipcsmF1 = 4       ! Fish consumption on species 1
      integer, parameter :: ipcsmF2 = 5       ! Fish consumption on species 2
      integer, parameter :: ipcsmF3 = 6       ! Fish consumption on species 3
      integer, parameter :: ipcsmF4 = 7       ! Fish consumption on species 4
      integer, parameter :: ipcsmF5 = 8       ! Fish consumption on species 5
! Species IDs
      integer, parameter :: ip_none = 0
      integer, parameter :: ip_albacore = 1
# endif
# ifdef FLEET
      integer, parameter :: ibiloc = 1         ! boat grid i location
      integer, parameter :: ibjloc = 2         ! boat grid j location
      integer, parameter :: ibport = 3         ! boat port
      integer, parameter :: ibnloc = 4         ! number of location visited
# endif
# ifdef EGGS_BISECTION
      logical :: lr_split = .true.
# endif
!
!  Other parameters.
!
      integer, dimension(Ngrids) :: idays
      integer, dimension(Ngrids) :: iyesterday
      integer, dimension(Ngrids) :: iyears
      integer, dimension(Ngrids) :: iyesteryear
!
!  Logical switches.
!
      logical, dimension(Ngrids) :: Lfish
# ifdef PREDATOR
      logical, dimension(Ngrids) :: Lpred
# endif
# ifdef FLEET
      logical, dimension(Ngrids) :: Lfleet
# endif

!
!  Fish restart switch.
!
      integer, dimension(Ngrids) :: fishrrec
# ifdef PREDATOR
      integer, dimension(Ngrids) :: predrrec
# endif
# ifdef FLEET
      integer, dimension(Ngrids) :: fleetrrec
# endif


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
      USE mod_biology
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

      allocate ( FISHES(ng) % bioenergy(NFishV(ng),Nfish(ng)) )

      allocate ( FISHES(ng) % feedback(NT(ng),Nfish(ng)) )

      allocate ( FISHES(ng) % species(Nfish(ng)) )

      allocate ( FISHES(ng) % lifestage(Nfish(ng)) )

      allocate ( FISHES(ng) % swimtype(Nfish(ng)) )

      allocate ( FISHES(ng) % alive(Nfish(ng)) )

      allocate ( FISHES(ng) % fishgirl(Nfish(ng)) )

      allocate ( FISHES(ng) % fishnodes(Nfish(ng)) )

      allocate ( FISHES(ng) % egg_dur(Nfish(ng)) )

      allocate ( FISHES(ng) % egg_num(Nfish(ng)) )

      allocate ( FISHES(ng) % ysac_dur(Nfish(ng)) )

      allocate ( FISHES(ng) % ysac_num(Nfish(ng)) )

      allocate ( FISHES(ng) % larv_dur(Nfish(ng)) )

      allocate ( FISHES(ng) % larv_num(Nfish(ng)) )

      allocate ( FISHES(ng) % juv_dur(Nfish(ng)) )

      allocate ( FISHES(ng) % juv_num(Nfish(ng)) )

      allocate ( FISHES(ng) % suba_num(Nfish(ng)) )

      allocate ( FISHES(ng) % sp_start(Nspecies(ng)) )

      allocate ( FISHES(ng) % sp_end(Nspecies(ng)) )

      allocate ( FISHES(ng) % fmortN(Nfish(ng)) )

      allocate ( FISHES(ng) % fmortS(Nfish(ng)) )

      allocate ( FISHES(ng) % fmortF(Nfish(ng)) )

      allocate ( FISHES(ng) % fmortP(Nfish(ng)) )

      allocate ( FISHES(ng) % cellid(Nfish(ng)) )

      allocate ( FISHES(ng) % num_free(Nspecies(ng)) )

      allocate ( FISHES(ng) % num_super(Nspecies(ng)) )

      allocate ( FISHES(ng) % next_free(Nspecies(ng)) )

      allocate ( FISHES(ng) % age_base(Nspecies(ng)) ) 

      allocate ( FISHES(ng) % species_base(Nspecies(ng)) )
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
        FISHES(ng) % fishgirl(iflt) = .FALSE.
        FISHES(ng) % fishnodes(iflt) % fish = iflt
        FISHES(ng) % fmortN(iflt) = IniVal
        FISHES(ng) % fmortS(iflt) = IniVal
        FISHES(ng) % fmortF(iflt) = IniVal
        FISHES(ng) % fmortP(iflt) = IniVal
        FISHES(ng) % cellid(iflt) = 0
      END DO

      Lfish(ng)=.TRUE.
# ifdef PREDATOR
      Lpred(ng)=.TRUE.
# endif
# ifdef FLEET
      Lfleet(ng)=.TRUE.
# endif
      LdefFISH(ng)=.TRUE.

      RETURN
      END SUBROUTINE allocate_fish
#endif
      END MODULE mod_fish
