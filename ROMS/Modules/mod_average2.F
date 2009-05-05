#include "cppdefs.h"
      MODULE mod_average2
#ifdef AVERAGES2
!
!svn $Id: mod_average.F 707 2008-08-19 17:58:55Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
! FOR the AVERAGE2 option, these become surface only fields.           !
!                                                                      !
!  2D Time-averaged fields for output purposes.                        !
!                                                                      !
!  avgu2d     2D velocity component (m/s) in the XI-direction.        !
!  avgv2d     2D velocity component (m/s) in the ETA-direction.       !
!  avgzeta    Free surface (m).                                       !
!                                                                      !
!  3D Time-averaged fields for output purposes.                        !
!                                                                      !
!  avglhf     Latent heat flux (W/m2).                                !
!  avglrf     Longwave radiation flux (W/m2).                         !
!  avgHuon    U-momentum flux, Hz*u/pn (m3/s).                        !
!  avgHuonT   Tracer u-transport, Hz*u*t/pn (Tunits m3/s).            !
!  avgHvom    V-momentum flux, Hz*v/pm (m3/s).                        !
!  avgHvomT   Tracer v-transport, Hz*v*t/pn (Tunits m3/s).            !
!  avgbus     Bottom u-momentum stress (N/m2).                        !
!  avgbvs     Bottom v-momentum stress (N/m2).                        !
!  avghbbl    Depth of oceanic bottom boundary layer (m).             !
!  avghsbl    Depth of oceanic surface boundary layer (m).            !
!  avgrho     Density anomaly (kg/m3).                                !
!  avgsssflx  Sea surface salinity flux correction.                   !
!  avgshf     Sensible heat flux (W/m2).                              !
!  avgsrf     Shortwave radiation flux (W/m2).                        !
!  avgstf     Surface net heat flux (W/m2).                           !
!  avgswf     Surface net salt flux (kg/m2/s).                        !
!  avgevap    Surface net evaporation (kg/m2/s).                      !
!  avgrain    Surface net rain fall (kg/m2/s).                        !
!  avgsus     Surface u-momentum stress (N/m2).                       !
!  avgsvs     Surface v-momentum stress (N/m2).                       !
!  avgt       Tracer type variables (usually, potential temperature   !
!               and salinity).                                         !
!  avguwind   2D wind velocity component (m/s) in the XI-direction.   !
!  avgvwind   2D wind velocity component (m/s) in the ETA-direction.  !
!  avgu3d     3D velocity component (m/s) in the XI-direction.        !
!  avgv3d     3D velocity component (m/s) in the ETA-direction.       !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_AVERAGE2

          real(r8), pointer :: avgu2d(:,:)
          real(r8), pointer :: avgv2d(:,:)
          real(r8), pointer :: avgzeta(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: avgrho(:,:)
          real(r8), pointer :: avgt(:,:,:)
#  ifdef BIO_GOANPZ
          real(r8), pointer :: avgst(:,:,:)
#  endif
          real(r8), pointer :: avgu3d(:,:)
          real(r8), pointer :: avgv3d(:,:)
#  ifdef AVERAGES_FLUXES
          real(r8), pointer :: avgstf(:,:)
          real(r8), pointer :: avgswf(:,:)
#   ifdef BULK_FLUXES
          real(r8), pointer :: avglhf(:,:)
          real(r8), pointer :: avglrf(:,:)
          real(r8), pointer :: avgshf(:,:)
          real(r8), pointer :: avguwind(:,:)
          real(r8), pointer :: avgvwind(:,:)
#    ifdef EMINUSP
          real(r8), pointer :: avgevap(:,:)
          real(r8), pointer :: avgrain(:,:)
#    endif
#   endif
#   ifdef SHORTWAVE
          real(r8), pointer :: avgsrf(:,:)
#   endif
#  endif
#  ifdef LMD_BKPP
          real(r8), pointer :: avghbbl(:,:)
#  endif
#  ifdef LMD_SKPP
          real(r8), pointer :: avghsbl(:,:)
#  endif
#  ifdef WRT_SSSFLX
          real(r8), pointer :: avgsssflx(:,:)
#  endif
#  ifdef ICE_MODEL
          real(r8), pointer :: avguice(:,:)
          real(r8), pointer :: avgvice(:,:)
          real(r8), pointer :: avgaice(:,:)
          real(r8), pointer :: avghice(:,:)
          real(r8), pointer :: avgtice(:,:)
          real(r8), pointer :: avgtimid(:,:)
          real(r8), pointer :: avghsno(:,:)
          real(r8), pointer :: avgsfwat(:,:)
          real(r8), pointer :: avgiomflx(:,:)
          real(r8), pointer :: avgageice(:,:)
          real(r8), pointer :: avgsig11(:,:)
          real(r8), pointer :: avgsig12(:,:)
          real(r8), pointer :: avgsig22(:,:)
          real(r8), pointer :: avgT0mk(:,:)
          real(r8), pointer :: avgS0mk(:,:)
          real(r8), pointer :: avgchu_iw(:,:)
          real(r8), pointer :: avgutau_iw(:,:)
#  endif
# endif
# ifdef AVERAGES_FLUXES
          real(r8), pointer :: avgsus(:,:)
          real(r8), pointer :: avgsvs(:,:)
          real(r8), pointer :: avgbus(:,:)
          real(r8), pointer :: avgbvs(:,:)
# endif

        END TYPE T_AVERAGE2

        TYPE (T_AVERAGE2), allocatable :: AVERAGE2(:)

      CONTAINS

      SUBROUTINE allocate_average2 (ng, LBi, UBi, LBj, UBj)
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
      IF (ng.eq.1 ) allocate ( AVERAGE2(Ngrids) )
!
      allocate ( AVERAGE2(ng) % avgu2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgv2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgzeta(LBi:UBi,LBj:UBj) )
# ifdef AVERAGES_FLUXES
      allocate ( AVERAGE2(ng) % avgsus(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgsvs(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgbus(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgbvs(LBi:UBi,LBj:UBj) )
# endif
# ifdef SOLVE3D
      allocate ( AVERAGE2(ng) % avgrho(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgt(LBi:UBi,LBj:UBj,NT(ng)) )
#  ifdef BIO_GOANPZ
      allocate ( AVERAGE2(ng) % avgst(LBi:UBi,LBj:UBj,NTS(ng)) )
#  endif
      allocate ( AVERAGE2(ng) % avgu3d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgv3d(LBi:UBi,LBj:UBj) )
#  ifdef AVERAGES_FLUXES
      allocate ( AVERAGE2(ng) % avgstf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgswf(LBi:UBi,LBj:UBj) )
#   ifdef BULK_FLUXES
      allocate ( AVERAGE2(ng) % avglhf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avglrf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgshf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avguwind(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgvwind(LBi:UBi,LBj:UBj) )
#    ifdef EMINUSP
      allocate ( AVERAGE2(ng) % avgevap(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgrain(LBi:UBi,LBj:UBj) )
#    endif
#   endif
#   ifdef SHORTWAVE
      allocate ( AVERAGE2(ng) % avgsrf(LBi:UBi,LBj:UBj) )
#   endif
#  endif

#  ifdef LMD_BKPP
      allocate ( AVERAGE2(ng) % avghbbl(LBi:UBi,LBj:UBj) )
#  endif
#  ifdef LMD_SKPP
      allocate ( AVERAGE2(ng) % avghsbl(LBi:UBi,LBj:UBj) )
#  endif
#  ifdef WRT_SSSFLX
      allocate ( AVERAGE2(ng) % avgsssflx(LBi:UBi,LBj:UBj) )
#  endif
#  ifdef ICE_MODEL
      allocate ( AVERAGE2(ng) % avguice(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgvice(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgaice(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avghice(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgtice(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgtimid(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avghsno(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgsfwat(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgiomflx(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgageice(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgsig11(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgsig12(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgsig22(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgT0mk(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgS0mk(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgutau_iw(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgchu_iw(LBi:UBi,LBj:UBj) )
#  endif
# endif

      RETURN
      END SUBROUTINE allocate_average2

      SUBROUTINE initialize_average2 (ng, tile)
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
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
# ifdef SOLVE3D
      integer :: itrc, k
# endif

      real(r8), parameter :: IniVal = 0.0_r8

# include "set_bounds.h"
!
!  Set array initialization range.
!
# ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
# else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          AVERAGE2(ng) % avgu2d(i,j) = IniVal
          AVERAGE2(ng) % avgv2d(i,j) = IniVal
          AVERAGE2(ng) % avgzeta(i,j) = IniVal
# ifdef AVERAGES_FLUXES
          AVERAGE2(ng) % avgsus(i,j) = IniVal
          AVERAGE2(ng) % avgsvs(i,j) = IniVal
          AVERAGE2(ng) % avgbus(i,j) = IniVal
          AVERAGE2(ng) % avgbvs(i,j) = IniVal
# endif
# ifdef SOLVE3D
#  ifdef AVERAGES_FLUXES
          AVERAGE2(ng) % avgstf(i,j) = IniVal
          AVERAGE2(ng) % avgswf(i,j) = IniVal
#   ifdef BULK_FLUXES
          AVERAGE2(ng) % avglhf(i,j) = IniVal
          AVERAGE2(ng) % avglrf(i,j) = IniVal
          AVERAGE2(ng) % avgshf(i,j) = IniVal
          AVERAGE2(ng) % avguwind(i,j) = IniVal
          AVERAGE2(ng) % avgvwind(i,j) = IniVal
#    ifdef EMINUSP
          AVERAGE2(ng) % avgevap(i,j) = IniVal
          AVERAGE2(ng) % avgrain(i,j) = IniVal
#    endif
#   endif
#   ifdef SHORTWAVE
          AVERAGE2(ng) % avgsrf(i,j) = IniVal
#   endif
#  endif
#  ifdef LMD_BKPP
          AVERAGE2(ng) % avghbbl(i,j) = IniVal
#  endif
#  ifdef LMD_SKPP
          AVERAGE2(ng) % avghsbl(i,j) = IniVal
#  endif
#  ifdef WRT_SSSFLX
          AVERAGE2(ng) % avgsssflx(i,j) = IniVal
#  endif
# endif
        END DO
# ifdef SOLVE3D
        DO i=Imin,Imax
            AVERAGE2(ng) % avgrho(i,j) = IniVal
            AVERAGE2(ng) % avgu3d(i,j) = IniVal
            AVERAGE2(ng) % avgv3d(i,j) = IniVal
        END DO
        DO itrc=1,NT(ng)
          DO i=Imin,Imax
            AVERAGE2(ng) % avgt(i,j,itrc) = IniVal
          END DO
        END DO
#  ifdef BIO_GOANPZ
        DO itrc=1,NTS(ng)
          DO i=Imin,Imax
            AVERAGE2(ng) % avgst(i,j,itrc) = IniVal
          END DO
        END DO
#  endif
        DO i=Imin,Imax
#  ifdef ICE_MODEL
          AVERAGE2(ng) % avguice(i,j) = IniVal
          AVERAGE2(ng) % avgvice(i,j) = IniVal
          AVERAGE2(ng) % avgaice(i,j) = IniVal
          AVERAGE2(ng) % avghice(i,j) = IniVal
          AVERAGE2(ng) % avgtice(i,j) = IniVal
          AVERAGE2(ng) % avgtimid(i,j) = IniVal
          AVERAGE2(ng) % avghsno(i,j) = IniVal
          AVERAGE2(ng) % avgsfwat(i,j) = IniVal
          AVERAGE2(ng) % avgiomflx(i,j) = IniVal
          AVERAGE2(ng) % avgageice(i,j) = IniVal
          AVERAGE2(ng) % avgsig11(i,j) = IniVal
          AVERAGE2(ng) % avgsig12(i,j) = IniVal
          AVERAGE2(ng) % avgsig22(i,j) = IniVal
          AVERAGE2(ng) % avgT0mk(i,j) = IniVal
          AVERAGE2(ng) % avgS0mk(i,j) = IniVal
          AVERAGE2(ng) % avgutau_iw(i,j) = IniVal
          AVERAGE2(ng) % avgchu_iw(i,j) = IniVal
#  endif
# endif

        END DO
      END DO

      RETURN
      END SUBROUTINE initialize_average2
#endif
      END MODULE mod_average2
