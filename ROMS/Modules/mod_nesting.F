      MODULE mod_nesting
#ifdef NESTING
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module defines generic structures used for nesting, composed,  !
!  and mosaic grids.                                                   !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!-----------------------------------------------------------------------
!  Nesting identification index of variables to process.
!-----------------------------------------------------------------------
!
!  The following identification indices are used in "initial" or
!  "main2d/main3d" to specify the variables that are processed in
!  each sub-timestep section.  Negative indices are used in grid
!  refinement whereas positive indices are used in composite grids.
!
      integer, parameter :: ndxdy = -4    ! extract on_u and om_v
      integer, parameter :: ngetD = -3    ! extract donor grid data
      integer, parameter :: nputD = -2    ! fill contact points
      integer, parameter :: n2way = -1    ! fine to course coupling
!
      integer, parameter :: nFSIC =  1    ! free surface initialization
      integer, parameter :: n2dIC =  2    ! 2D momentum initialization
      integer, parameter :: n3dIC =  3    ! 3D momentum initialization
      integer, parameter :: nTVIC =  4    ! tracers initialization
      integer, parameter :: nbstr =  5    ! bottom stress (bustr,bvstr)
      integer, parameter :: nrhst =  6    ! RHS terms (tracers)
      integer, parameter :: nzeta =  7    ! 3D kernel free-surface
      integer, parameter :: nzwgt =  8    ! 3D vertical weights
      integer, parameter :: n2dPS =  9    ! 2D engine Predictor Step
      integer, parameter :: n2dCS = 10    ! 2D engine Corrector Step
      integer, parameter :: n2dfx = 11    ! time-averaged 2D fluxes
      integer, parameter :: n3duv = 12    ! 3D momentum and fluxes
      integer, parameter :: n3dTV = 13    ! 3D tracer variables
!
!-----------------------------------------------------------------------
!  Nesting parameters.
!-----------------------------------------------------------------------
!
!  Nested grid connectivity switches. It is used to determine the
!  dimensions of the numerical kernel allocatable arrays. The arrays
!  have extra points due to the contact regions in any of sides
!  of the physical grid (1=iwest, 2=isouth, 3=ieast, 4=inorth).
!
      logical, allocatable :: ContactRegion(:,:)    ! [4,Ngrids]
!
!  Data donor REFINED structure index associated with refinement grid.
!  Recall that both refinement and composite grids are allow and this
!  mapping is for efficient memory usage.
!
      integer, allocatable :: RefineDonor(:)        ! [Ngrids]
!
!  Donor and reciver grids for each contact region. These paremeters
!  are aslo duplicated in the T_NGC structure.
!
      integer, allocatable :: donor_grid(:)         ! [Ncontact]
      integer, allocatable :: receiver_grid(:)      ! [Ncontact]
!
!  Compact arrays used to unpack data from nested grids contact points
!  NetCDF file. They are allocated to the size "datum" dimension in
!  routine "set_contact".  The start and end indices for each C-type
!  variable are used to unpack from compact vector.
!
      integer :: NCdatum
      integer, allocatable :: NCpoints(:)           ! [Ncontact]
      integer, allocatable :: NstrR(:), NendR(:)    ! [Ncontact]
      integer, allocatable :: NstrU(:), NendU(:)    ! [Ncontact]
      integer, allocatable :: NstrV(:), NendV(:)    ! [Ncontact]
!
      integer, allocatable :: contact_region(:)     ! [NCdatum]
      integer, allocatable :: on_boundary(:)        ! [NCdatum]
      integer, allocatable :: Idg_cp(:)             ! [NCdatum]
      integer, allocatable :: Jdg_cp(:)             ! [NCdatum]
      integer, allocatable :: Irg_cp(:)             ! [NCdatum]
      integer, allocatable :: Jrg_cp(:)             ! [NCdatum]
!
! Generic structure for 2D state fields.
!
        TYPE T_NEST2D
          integer, pointer :: id             ! state field id
          integer, pointer :: gtype          ! staggered grid type
          integer, pointer :: tile(:)        ! tile partition

          integer, pointer :: LBi(:)         ! lower bound I-dimension
          integer, pointer :: UBi(:)         ! upper bound I-dimension
          integer, pointer :: LBj(:)         ! lower bound J-dimension
          integer, pointer :: UBj(:)         ! upper bound J-dimension

          integer, pointer :: Istr(:)        ! starting tile I-direction
          integer, pointer :: Iend(:)        ! ending   tile I-direction
          integer, pointer :: Jstr(:)        ! starting tile J-direction
          integer, pointer :: Jend(:)        ! ending   tile J-direction

          integer, pointer :: IstrR(:)       ! starting RHO  I-direction
          integer, pointer :: IstrU(:)       ! starting U    I-direction
          integer, pointer :: IendR(:)       ! ending   RHO  I-direction

          integer, pointer :: JstrR(:)       ! starting RHO  J-direction
          integer, pointer :: JstrV(:)       ! starting V    J-direction
          integer, pointer :: JendR(:)       ! ending   RHO  J-direction

          real (r8), pointer :: x(:,:)       ! X-positions
          real (r8), pointer :: y(:,:)       ! Y-positions
# ifdef MASKING
          real (r8), pointer :: mask(:,:)    ! land-sea masking
# endif
          real (r8), pointer :: s(:,:)       ! state array(i,j)
        END TYPE T_NEST2D

# ifdef SOLVE3D
!
! Generic structure for 3D state fields.
!
        TYPE T_NEST3D
          integer, pointer :: id             ! state field id
          integer, pointer :: gtype          ! staggered grid type
          integer, pointer :: tile(:)        ! tile partition

          integer, pointer :: LBi(:)         ! lower bound I-dimension
          integer, pointer :: UBi(:)         ! upper bound I-dimension
          integer, pointer :: LBj(:)         ! lower bound J-dimension
          integer, pointer :: UBj(:)         ! upper bound J-dimension
          integer, pointer :: LBk            ! lower bound K-dimension
          integer, pointer :: UBk            ! upper bound K-dimension

          integer, pointer :: Istr(:)        ! starting tile I-direction
          integer, pointer :: Iend(:)        ! ending   tile I-direction
          integer, pointer :: Jstr(:)        ! starting tile J-direction
          integer, pointer :: Jend(:)        ! ending   tile J-direction

          integer, pointer :: IstrR(:)       ! starting RHO  I-direction
          integer, pointer :: IstrU(:)       ! starting U    I-direction
          integer, pointer :: IendV(:)       ! ending   RHO  I-direction

          integer, pointer :: JstrR(:)       ! starting RHO  J-direction
          integer, pointer :: JstrV(:)       ! starting V    J-direction
          integer, pointer :: JendR(:)       ! ending   RHO  J-direction

          integer, pointer :: Kstr           ! starting K-index
          integer, pointer :: Kend           ! ending   K-index

          real (r8), pointer :: x(:,:)       ! X-positions
          real (r8), pointer :: y(:,:)       ! Y-positions
#  ifdef MASKING
          real (r8), pointer :: mask(:,:)    ! land-sea masking
#  endif
          real (r8), pointer :: s(:,:,:)     ! state array(i,j,k)
        END TYPE T_NEST3D
!
! Generic structure for 4D state fields.
!
        TYPE T_NEST4D
          integer, pointer :: id             ! state field id
          integer, pointer :: gtype          ! staggered grid type
          integer, pointer :: tile(:)        ! tile partition

          integer, pointer :: LBi(:)         ! lower bound I-dimension
          integer, pointer :: UBi(:)         ! upper bound I-dimension
          integer, pointer :: LBj(:)         ! lower bound J-dimension
          integer, pointer :: UBj(:)         ! upper bound J-dimension
          integer, pointer :: LBk            ! lower bound K-dimension
          integer, pointer :: UBk            ! upper bound K-dimension
          integer, pointer :: LBl            ! lower bound L-dimension
          integer, pointer :: UBl            ! upper bound L-dimension

          integer, pointer :: Istr(:)        ! starting tile I-direction
          integer, pointer :: Iend(:)        ! ending   tile I-direction
          integer, pointer :: Jstr(:)        ! starting tile J-direction
          integer, pointer :: Jend(:)        ! ending   tile J-direction

          integer, pointer :: IstrR(:)       ! starting RHO  I-direction
          integer, pointer :: IstrU(:)       ! starting U    I-direction
          integer, pointer :: IendV(:)       ! ending   RHO  I-direction

          integer, pointer :: JstrR(:)       ! starting RHO  J-direction
          integer, pointer :: JstrV(:)       ! starting V    J-direction
          integer, pointer :: JendR(:)       ! ending   RHO  J-direction

          integer, pointer :: Kstr           ! starting K-index
          integer, pointer :: Kend           ! ending   K-index
          integer, pointer :: Lstr           ! starting L-index
          integer, pointer :: Lend           ! ending   L-index

          real (r8), pointer :: x(:,:)       ! X-positions
          real (r8), pointer :: y(:,:)       ! Y-positions
#  ifdef MASKING
          real (r8), pointer :: mask(:,:)    ! land-sea masking
#  endif
          real (r8), pointer :: s(:,:,:,:)   ! state array(i,j,k,l)
        END TYPE T_NEST4D
# endif

      CONTAINS

      SUBROUTINE allocate_nesting2d (field2d, id, gtype, Tindex)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 2D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, target :: id, gtype

      integer, intent(in), optional :: Tindex

      TYPE (T_NEST2D), pointer :: field2d(:)
!
!  Local variable declarations.
!
      integer :: ng
!
!-----------------------------------------------------------------------
!  Generic 2D structure.
!-----------------------------------------------------------------------
!
!  Allocate.
!
      IF (ASSOCIATED(field2d)) THEN
        deallocate ( field2d )
      END IF
      allocate ( field2d(Ngrids) )
!
!  Load field descriptors.
!
      DO ng=1,Ngrids
!
!  Load id, grid type, tile, array bounds, and starting/ending
!  computational indices.
!
        field2d(ng) % id    => id
        field2d(ng) % gtype => gtype
        field2d(ng) % tile  => BOUNDS(ng) % tile
        field2d(ng) % LBi   => BOUNDS(ng) % LBi
        field2d(ng) % UBi   => BOUNDS(ng) % UBi
        field2d(ng) % LBj   => BOUNDS(ng) % LBj
        field2d(ng) % UBj   => BOUNDS(ng) % UBj
        field2d(ng) % Istr  => BOUNDS(ng) % Istr
        field2d(ng) % Iend  => BOUNDS(ng) % Iend
        field2d(ng) % Jstr  => BOUNDS(ng) % Jstr
        field2d(ng) % Jend  => BOUNDS(ng) % Jend
        field2d(ng) % IstrR => BOUNDS(ng) % IstrR
        field2d(ng) % IstrU => BOUNDS(ng) % IstrU
        field2d(ng) % IendR => BOUNDS(ng) % IendR
        field2d(ng) % JstrR => BOUNDS(ng) % JstrR
        field2d(ng) % JstrV => BOUNDS(ng) % JstrV
        field2d(ng) % JendR => BOUNDS(ng) % JendR
!
!  Associate the appropriate grid arrays.
!
        IF (spherical) THEN
          IF (gtype.eq.u2dvar) THEN
            field2d(ng) % x => GRID(ng) % lonu
            field2d(ng) % y => GRID(ng) % latu
          ELSE IF (gtype.eq.v2dvar) THEN
            field2d(ng) % x => GRID(ng) % lonv
            field2d(ng) % y => GRID(ng) % latv
          ELSE
            field2d(ng) % x => GRID(ng) % lonr
            field2d(ng) % y => GRID(ng) % latr
          END IF
        ELSE
          IF (gtype.eq.u2dvar) THEN
            field2d(ng) % x => GRID(ng) % xu
            field2d(ng) % y => GRID(ng) % yu
          ELSE IF (gtype.eq.v2dvar) THEN
            field2d(ng) % x => GRID(ng) % xv
            field2d(ng) % y => GRID(ng) % yv
          ELSE
            field2d(ng) % x => GRID(ng) % xr
            field2d(ng) % y => GRID(ng) % yr
          END IF
        END IF

# ifdef MASKING
!
!  Associate the appropriate Land/Sea mask.
!
        IF (gtype.eq.u2dvar) THEN
          field2d(ng) % mask  => GRID(ng) % umask
        ELSE IF (gtype.eq.v2dvar) THEN
          field2d(ng) % mask  => GRID(ng) % vmask
        ELSE
          field2d(ng) % mask  => GRID(ng) % rmask
        END IF
# endif
!
!  Associate the appropriate state 2D array.
!
        IF (id.eq.idFsur) THEN
          field2d(ng) % s => OCEAN(ng) % zeta(:,:,Tindex)
        ELSE IF (id.eq.idUbar) THEN
          field2d(ng) % s => OCEAN(ng) % ubar(:,:,Tindex)
        ELSE IF (id.eq.idVbar) THEN
          field2d(ng) % s => OCEAN(ng) % vbar(:,:,Tindex)
        ELSE IF (id.eq.idRzet) THEN
          field2d(ng) % s => OCEAN(ng) % rzeta(:,:,Tindex)
# ifdef SOLVE3D
        ELSE IF (id.eq.idZavg) THEN
          field2d(ng) % s => COUPLING(ng) % Zt_avg1
        ELSE IF (id.eq.idUfx1) THEN
          field2d(ng) % s => COUPLING(ng) % DU_avg1
        ELSE IF (id.eq.idVfx1) THEN
          field2d(ng) % s => COUPLING(ng) % DV_avg1
# endif
        END IF

      END DO

      RETURN
      END SUBROUTINE allocate_nesting2d

        real(r8), pointer :: on_u(:)          ! [Npoints]
        real(r8), pointer :: om_v(:)          ! [Npoints]

# ifdef SOLVE3D

      SUBROUTINE allocate_nesting3d (field3d, id, gtype, Tindex)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 3D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Imin, Imax, Jmin, Jmax
      integer :: cr, dg, ng, rg
      integer :: i, ibry, ic, j, m, my_tile
      integer :: ispval

      TYPE (T_NEST3D), pointer :: field3d(:)
!
!  Local variable declarations.
!
!  Allocate boundary connectivity (type T_BCP) structure.
!
      allocate ( BRY_CONTACT(4,Ncontact) )
!
!  Allocate arrays in boundary connectivity structure.
!
      my_tile=-1                           ! for global values
      DO cr=1,Ncontact
        dg=receiver_grid(cr)
        LBi=BOUNDS(dg)%LBi(my_tile)
        UBi=BOUNDS(dg)%UBi(my_tile)
        LBj=BOUNDS(dg)%LBj(my_tile)
        UBj=BOUNDS(dg)%UBj(my_tile)
        DO ibry=1,4
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              allocate ( BRY_CONTACT(ibry,cr) % Ib(LBj:UBj) )
              allocate ( BRY_CONTACT(ibry,cr) % Jb(LBj:UBj) )
              allocate ( BRY_CONTACT(ibry,cr) % C2Bindex(LBj:UBj) )
              allocate ( BRY_CONTACT(ibry,cr) % Tflux(LBj:UBj,          &
     &                                                N(dg),NT(dg)) )
            CASE (isouth, inorth)
              allocate ( BRY_CONTACT(ibry,cr) % Ib(LBi:UBi) )
              allocate ( BRY_CONTACT(ibry,cr) % Jb(LBi:UBi) )
              allocate ( BRY_CONTACT(ibry,cr) % C2Bindex(LBi:UBi) )
              allocate ( BRY_CONTACT(ibry,cr) % Tflux(LBi:UBi,          &
     &                                                N(dg),NT(dg)) )
          END SELECT
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Generic 3D structure.
!-----------------------------------------------------------------------
!
!  Allocate.
!
      IF (ASSOCIATED(field3d)) THEN
        deallocate ( field3d )
      END IF
      IF (.not.allocated(Jbmax)) THEN
        allocate ( Jbmax(4,Ncontact) )
        Jbmax = ispval
      END IF
!
      DO cr=1,Ncontact
        DO ibry=1,4
          BRY_CONTACT(ibry,cr) % spv = ispval
          BRY_CONTACT(ibry,cr) % Ib = ispval
          BRY_CONTACT(ibry,cr) % Jb = ispval
          BRY_CONTACT(ibry,cr) % C2Bindex = ispval
        END DO
      END DO
!
!  Identify contact points located on the grid boundary.  Notice that
!  all the contact points are looped "Ncontact" times to determine the
!  cell containing the finer grid physical boundary for both donor and
!  receiver grids. In refinement grids, the contact regions are flipped
!  for consistency with the coding logic elsewhere.
!
      DO cr=1,Ncontact
        dg=donor_grid(cr)
        rg=receiver_grid(cr)
        ic=1
        DO m=1,NCdatum
          ibry=on_boundary(m)
          IF ((m.gt.NCpoints(ic)).and.(ic.lt.Ncontact)) THEN
            ic=ic+1
          END IF
          IF (RefinedGrid(dg).and.RefinedGrid(rg)) THEN
            IF (rg.lt.dg) THEN
              i=Idg_cp(m)        ! coarser donor grid cell containing
              j=Jdg_cp(m)        ! finer grid physical boundary points
            ELSE
              i=Irg_cp(m)        ! finer receiver grid cell physical
              j=Jrg_cp(m)        ! boundary points
            END IF
          ELSE
            i=Irg_cp(m)          ! composite grids
            j=Jrg_cp(m)
          END IF
          IF ((ibry.eq.iwest ).or.(ibry.eq.ieast )) THEN
            Ibmin(ibry,cr)=MIN(i,Ibmin(ibry,cr))
            Ibmax(ibry,cr)=MAX(i,Ibmax(ibry,cr))
            Jbmin(ibry,cr)=MIN(j,Jbmin(ibry,cr))
            Jbmax(ibry,cr)=MAX(j,Jbmax(ibry,cr))
            BRY_CONTACT(ibry,cr) % Ib(j) = i
            BRY_CONTACT(ibry,cr) % Jb(j) = j
            BRY_CONTACT(ibry,cr) % C2Bindex(j) = m-NstrU(ic)+1
          ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
            Ibmin(ibry,cr)=MIN(i,Ibmin(ibry,cr))
            Ibmax(ibry,cr)=MAX(i,Ibmax(ibry,cr))
            Jbmin(ibry,cr)=MIN(j,Jbmin(ibry,cr))
            Jbmax(ibry,cr)=MAX(j,Jbmax(ibry,cr))
            BRY_CONTACT(ibry,cr) % Ib(i) = i
            BRY_CONTACT(ibry,cr) % Jb(i) = j
            BRY_CONTACT(ibry,cr) % C2Bindex(i) = m-NstrV(ic)+1
          END IF
        END DO
      END DO
!
!  Associate the appropriate grid arrays.
!
        IF (spherical) THEN
          IF (gtype.eq.u3dvar) THEN
            field3d(ng) % x => GRID(ng) % lonu
            field3d(ng) % y => GRID(ng) % latu
          ELSE IF (gtype.eq.v3dvar) THEN
            field3d(ng) % x => GRID(ng) % lonv
            field3d(ng) % y => GRID(ng) % latv
          ELSE
            field3d(ng) % x => GRID(ng) % lonr
            field3d(ng) % y => GRID(ng) % latr
          END IF
        ELSE
          IF (gtype.eq.u3dvar) THEN
            field3d(ng) % x => GRID(ng) % xu
            field3d(ng) % y => GRID(ng) % yu
          ELSE IF (gtype.eq.v3dvar) THEN
            field3d(ng) % x => GRID(ng) % xv
            field3d(ng) % y => GRID(ng) % yv
          ELSE
            field3d(ng) % x => GRID(ng) % xr
            field3d(ng) % y => GRID(ng) % yr
          END IF
        END IF

#  ifdef MASKING
!
!  Associate the appropriate Land/Sea mask.
!
        IF (gtype.eq.u3dvar) THEN
          field3d(ng) % mask  => GRID(ng) % umask
        ELSE IF (gtype.eq.v3dvar) THEN
          field3d(ng) % mask  => GRID(ng) % vmask
        ELSE
          field3d(ng) % mask  => GRID(ng) % rmask
        END IF
#  endif
!
!  Associate the appropriate state 3D array.
!
        IF (id.eq.idUvel) THEN
          field3d(ng) % s => OCEAN(ng) % u(:,:,:,Tindex)
        ELSE IF (id.eq.idVvel) THEN
          field3d(ng) % s => OCEAN(ng) % v(:,:,:,Tindex)
        ELSE IF (id.eq.idRu3d) THEN
          field3d(ng) % s => OCEAN(ng) % ru(:,:,:,Tindex)
        ELSE IF (id.eq.idRv3d) THEN
          field3d(ng) % s => OCEAN(ng) % rv(:,:,:,Tindex)
        ELSE IF (id.eq.idOvel) THEN
          field3d(ng) % s => OCEAN(ng) % W
        ELSE IF (id.eq.idWvel) THEN
          field3d(ng) % s => OCEAN(ng) % wvel
        ELSE IF (id.eq.idDano) THEN
          field3d(ng) % s => OCEAN(ng) % rho
        ELSE IF (id.eq.idVvis) THEN
          field3d(ng) % s => MIXING(ng) % Akv
#  if defined GLS_MIXING || defined MY25_MIXING
        ELSE IF (id.eq.idMtke) THEN
          field3d(ng) % s => MIXING(ng) % tke(:,:,:,Tindex)
        ELSE IF (id.eq.idMtls) THEN
          field3d(ng) % s => MIXING(ng) % gls(:,:,:,Tindex)
        ELSE IF (id.eq.idVmKK) THEN
          field3d(ng) % s => MIXING(ng) % Akk
#   ifdef GLS_MIXING
        ELSE IF (id.eq.idVvel) THEN
          field3d(ng) % s => MIXING(ng) % Akp
#   endif
#  endif
        END IF

      END DO

      RETURN
      END SUBROUTINE allocate_nesting3d

      SUBROUTINE allocate_nesting4d (field4d, id, gtype, Tindex)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 4D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, target :: id, gtype

      integer, intent(in), optional :: Tindex

      TYPE (T_NEST4D), pointer :: field4d(:)
!
!  Local variable declarations.
!
      integer :: ng
      integer, target :: LBk, UBk, LBl, UBl
      integer, target :: Kstr, Kend, Lstr, Lend
!
!-----------------------------------------------------------------------
!  Generic 4D structure.
!-----------------------------------------------------------------------
!
!  Allocate.
!
      IF (associated(field4d)) THEN
        deallocate ( field4d )
      END IF
      allocate ( field4d(Ngrids) )
!
!  Load field descriptors.
!
      DO ng=1,Ngrids
!
!  Load id, grid type, tile, array bounds, and starting/ending
!  computational indices.
!
        field4d(ng) % id    => id
        field4d(ng) % gtype => gtype
        field4d(ng) % tile  => BOUNDS(ng) % tile
        field4d(ng) % LBi   => BOUNDS(ng) % LBi
        field4d(ng) % UBi   => BOUNDS(ng) % UBi
        field4d(ng) % LBj   => BOUNDS(ng) % LBj
        field4d(ng) % UBj   => BOUNDS(ng) % UBj
        field4d(ng) % Istr  => BOUNDS(ng) % Istr
        field4d(ng) % Iend  => BOUNDS(ng) % Iend
        field4d(ng) % Jstr  => BOUNDS(ng) % Jstr
        field4d(ng) % Jend  => BOUNDS(ng) % Jend
        field4d(ng) % IstrR => BOUNDS(ng) % IstrR
        field4d(ng) % IstrU => BOUNDS(ng) % IstrU
        field4d(ng) % IendR => BOUNDS(ng) % IendR
        field4d(ng) % JstrR => BOUNDS(ng) % JstrR
        field4d(ng) % JstrV => BOUNDS(ng) % JstrV
        field4d(ng) % JendR => BOUNDS(ng) % JendR
!
!  Load third-dimension bounds.
!
        IF (gtype.eq.r3dvar) THEN
          LBk=1
          UBk=N(ng)
          Kstr=1
          Kend=N(ng)
          field4d(ng) % LBk   => LBk
          field4d(ng) % UBk   => UBk
          field4d(ng) % Kstr  => Kstr
          field4d(ng) % Kend  => Kend
        ELSE IF
        END IF
!
!  Associate the appropriate grid arrays.
!
        IF (spherical) THEN
          IF (gtype.eq.u3dvar) THEN
            field4d(ng) % x => GRID(ng) % lonu
            field4d(ng) % y => GRID(ng) % latu
          ELSE IF (gtype.eq.v3dvar) THEN
            field4d(ng) % x => GRID(ng) % lonv
            field4d(ng) % y => GRID(ng) % latv
          ELSE
            field4d(ng) % x => GRID(ng) % lonr
            field4d(ng) % y => GRID(ng) % latr
          END IF
        ELSE
          IF (gtype.eq.u3dvar) THEN
            field4d(ng) % x => GRID(ng) % xu
            field4d(ng) % y => GRID(ng) % yu
          ELSE IF (gtype.eq.v3dvar) THEN
            field4d(ng) % x => GRID(ng) % xv
            field4d(ng) % y => GRID(ng) % yv
          ELSE
            field4d(ng) % x => GRID(ng) % xr
            field4d(ng) % y => GRID(ng) % yr
          END IF
        END IF

#  ifdef MASKING
!
!  Associate the appropriate Land/Sea mask.
!
        IF (gtype.eq.u3dvar) THEN
          field4d(ng) % mask  => GRID(ng) % umask
        ELSE IF (gtype.eq.v3dvar) THEN
          field4d(ng) % mask  => GRID(ng) % vmask
        ELSE
          field4d(ng) % mask  => GRID(ng) % rmask
        END IF
#  endif
!
!  Associate the appropriate state 4D array.  Notice that temperature
!  triggers associating all tracers.
!
        IF (id.eq.idTvar(itemp)) THEN
          LBl=1
          UBl=NT(ng)
          Lstr=1
          Lend=NT(ng)
          field4d(ng) % LBl   => LBl
          field4d(ng) % UBl   => UBl
          field4d(ng) % Lstr  => Lstr
          field4d(ng) % Lend  => Lend
          field4d(ng) % s     => OCEAN(ng) % t(:,:,:,Tindex,:)
        END IF

          REFINED(cr) % on_u = IniVal
          REFINED(cr) % om_v = IniVal

# ifdef SOLVE3D
          REFINED(cr) % u = IniVal
          REFINED(cr) % v = IniVal

      RETURN
      END SUBROUTINE allocate_nesting4d
# endif

#endif
      END MODULE mod_nesting
