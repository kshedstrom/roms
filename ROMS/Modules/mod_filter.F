#include "cppdefs.h"
      MODULE mod_filter
#ifdef FILTERED
!
!=======================================================================
!  Copyright (c) 2002-2009 ROMS/TOMS Group                             !
!================================================== Hernan G. Arango ===
!                                                                      !
!  2D Time-filtered fields for output purposes.                        !
!                                                                      !
!  filu2d     2D velocity component (m/s) in the XI-direction.         !
!  filv2d     2D velocity component (m/s) in the ETA-direction.        !
!  filzeta    Free surface (m).                                        !
!                                                                      !
!  3D Time-filtered fields for output purposes.                        !
!                                                                      !
!  filrho     Density anomaly (kg/m3).                                 !
!  filt       Tracer type variables (usually, potential temperature    !
!               and salinity).                                         !
!  filu3d     3D velocity component (m/s) in the XI-direction.         !
!  filv3d     3D velocity component (m/s) in the ETA-direction.        !
!  filw3d     S-coordinate [omega*Hz/mn] vertical velocity (m3/s).     !
!                                                                      !
!  2D Time-filtered ice fields for output purposes.                    !
!                                                                      !
!  filrho     Density anomaly (kg/m3).                                 !
!                                                                      !
!  Filter characteristics                                              !
!                                                                      !
!  FIL        Number of concurrent filter sums.                        !
!  dtfil      Filter sampling interval in seconds.                     !
!  nw         Number of filter weights.                                !
!  nfil       Filter length in timesteps.                              !
!  fcount     Tracks how many filtering intervals have been summed.    !
!  fw         Filter weights.                                          !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_FILTER

          real(r8), pointer :: filzeta(:,:,:)
          real(r8), pointer :: filu2d(:,:,:)
          real(r8), pointer :: filv2d(:,:,:)
# ifdef SOLVE3D
          real(r8), pointer :: filu3d(:,:,:,:)
          real(r8), pointer :: filv3d(:,:,:,:)
#  ifndef FILTRIM
          real(r8), pointer :: filrho(:,:,:,:)
          real(r8), pointer :: filt(:,:,:,:,:)
          real(r8), pointer :: filw3d(:,:,:,:)
#  endif
#  ifdef LMD_SKPP
          real(r8), pointer :: filhsbl(:,:,:)
#  endif
#  ifdef LMD_BKPP
          real(r8), pointer :: filhbbl(:,:,:)
#  endif
#  ifdef AVERAGES_AKV
          real(r8), pointer :: filAKv(:,:,:,:)
#  endif
#  ifdef AVERAGES_AKT
          real(r8), pointer :: filAKt(:,:,:,:)
#  endif
#  ifdef AVERAGES_AKS
          real(r8), pointer :: filAKs(:,:,:,:)
#  endif
#  ifdef AVERAGES_FLUXES
          real(r8), pointer :: filstf(:,:,:)
          real(r8), pointer :: filswf(:,:,:)
#   ifdef BULK_FLUXES
          real(r8), pointer :: fillhf(:,:,:)
          real(r8), pointer :: filshf(:,:,:)
          real(r8), pointer :: fillrf(:,:,:)
#   endif
#   ifdef SHORTWAVE
          real(r8), pointer :: filsrf(:,:,:)
#   endif
#  endif
#  ifdef ICE_MODEL
          real(r8), pointer :: filui(:,:,:)
          real(r8), pointer :: filvi(:,:,:)
          real(r8), pointer :: filai(:,:,:)
          real(r8), pointer :: filhi(:,:,:)
          real(r8), pointer :: filhsn(:,:,:)
          real(r8), pointer :: filsfwat(:,:,:)
          real(r8), pointer :: filtis(:,:,:)
          real(r8), pointer :: filti(:,:,:)
          real(r8), pointer :: filAgeice(:,:,:)
          real(r8), pointer :: fils0mk(:,:,:)
          real(r8), pointer :: filt0mk(:,:,:)
          real(r8), pointer :: filWfr(:,:,:)
          real(r8), pointer :: filWai(:,:,:)
          real(r8), pointer :: filWao(:,:,:)
          real(r8), pointer :: filWio(:,:,:)
          real(r8), pointer :: filWro(:,:,:)
          real(r8), pointer :: filsig11(:,:,:)
          real(r8), pointer :: filsig22(:,:,:)
          real(r8), pointer :: filsig12(:,:,:)
          real(r8), pointer :: filio_mflux(:,:,:)
          real(r8), pointer :: filutau_iw(:,:,:)
          real(r8), pointer :: filchu_iw(:,:,:)
#  endif
# endif
# ifdef AVERAGES_FLUXES
          real(r8), pointer :: filsus(:,:,:)
          real(r8), pointer :: filsvs(:,:,:)
# endif

        END TYPE T_FILTER

        TYPE (T_FILTER), allocatable :: FILTER(:)

!    Hourly filter weights for pl33 filter.
#undef PL33
!    Hourly filter weights for 5 day lolopass2 filter
#define LOLOPASS2
# ifdef PL33
        integer, parameter :: FIL = 3
        integer, parameter :: dtfil = 3600
        integer, parameter :: nw = 72
        integer :: nfil
        integer, dimension(FIL) :: fcount
        real(r8), dimension(nw) :: fw =                                 &
     &            (/                                                    &
     &             0.0,      0.0,      0.0,      0.0,      0.0,         &
     &            -0.00027, -0.00114, -0.00211, -0.00317, -0.00427,     &
     &            -0.00537, -0.00641, -0.00735, -0.00811, -0.00864,     &
     &            -0.00887, -0.00872, -0.00816, -0.00714, -0.00560,     &
     &            -0.00355, -0.00097,  0.00213,  0.00574,  0.00980,     &
     &             0.01425,  0.01902,  0.02400,  0.02911,  0.03423,     &
     &             0.03923,  0.04399,  0.04842,  0.05237,  0.05576,     &
     &             0.05850,  0.06051,  0.06174,  0.06215,               &
     &             0.06174,  0.06051,  0.05850,                         &
     &             0.05576,  0.05237,  0.04842,  0.04399,  0.03923,     &
     &             0.03423,  0.02911,  0.02400,  0.01902,  0.01425,     &
     &             0.00980,  0.00574,  0.00213, -0.00097, -0.00355,     &
     &            -0.00560, -0.00714, -0.00816, -0.00872, -0.00887,     &
     &            -0.00864, -0.00811, -0.00735, -0.00641, -0.00537,     &
     &            -0.00427, -0.00317, -0.00211, -0.00114, -0.00027      &
     &              /)
# elif defined LOLOPASS2
        integer, parameter :: FIL = 5
        integer, parameter :: dtfil = 3600
        integer, parameter :: nw = 120
        integer :: nfil
        integer, dimension(FIL) :: fcount
        real(r8), dimension(nw) :: fw =                                 &
     &   (/                                                             &
     &     0.00000e+0, 2.24482e-6, 1.01382e-5, 2.69017e-5, 5.52832e-5,  &
     &     9.72726e-5, 1.53843e-4, 2.24728e-4, 3.08265e-4, 4.01287e-4,  &
     &     4.99108e-4, 5.95586e-4, 6.83272e-4, 7.53646e-4, 7.97444e-4,  &
     &     8.05050e-4, 7.66951e-4, 6.74250e-4, 5.19196e-4, 2.95727e-4,  &
     &     0.00000e+0,-3.69122e-4,-8.09640e-4,-1.31607e-3,-1.87920e-3,  &
     &    -2.48599e-3,-3.11952e-3,-3.75924e-3,-4.38118e-3,-4.95844e-3,  &
     &    -5.46179e-3,-5.86035e-3,-6.12246e-3,-6.21655e-3,-6.11216e-3,  &
     &    -5.78090e-3,-5.19754e-3,-4.34092e-3,-3.19491e-3,-1.74920e-3,  &
     &     0.00000e+0, 2.04941e-3, 4.38837e-3, 6.99864e-3, 9.85452e-3,  &
     &     1.29230e-2, 1.61646e-2, 1.95334e-2, 2.29786e-2, 2.64454e-2,  &
     &     2.98760e-2, 3.32112e-2, 3.63914e-2, 3.93586e-2, 4.20574e-2,  &
     &     4.44366e-2, 4.64503e-2, 4.80592e-2, 4.92317e-2, 4.99447e-2,  &
     &     5.01840e-2,                                                  &
     &     4.99447e-2, 4.92317e-2, 4.80592e-2, 4.64503e-2, 4.44366e-2,  &
     &     4.20574e-2, 3.93586e-2, 3.63914e-2, 3.32112e-2, 2.98760e-2,  &
     &     2.64454e-2, 2.29786e-2, 1.95334e-2, 1.61646e-2, 1.29230e-2,  &
     &     9.85452e-3, 6.99864e-3, 4.38837e-3, 2.04941e-3, 0.00000e+0,  &
     &    -1.74920e-3,-3.19491e-3,-4.34092e-3,-5.19754e-3,-5.78090e-3,  &
     &    -6.11216e-3,-6.21655e-3,-6.12246e-3,-5.86035e-3,-5.46179e-3,  &
     &    -4.95844e-3,-4.38118e-3,-3.75924e-3,-3.11952e-3,-2.48599e-3,  &
     &    -1.87920e-3,-1.31607e-3,-8.09640e-4,-3.69122e-4, 0.00000e+0,  &
     &     2.95727e-4, 5.19196e-4, 6.74250e-4, 7.66951e-4, 8.05050e-4,  &
     &     7.97444e-4, 7.53646e-4, 6.83272e-4, 5.95586e-4, 4.99108e-4,  &
     &     4.01287e-4, 3.08265e-4, 2.24728e-4, 1.53843e-4, 9.72726e-5,  &
     &     5.52832e-5, 2.69017e-5, 1.01382e-5, 2.24482e-6               &
     &              /)
# else
        integer, parameter :: FIL = 2
        integer, parameter :: dtfil = 900
        integer, parameter :: nw = 4
        integer :: nfil
        integer, dimension(FIL) :: fcount
        real(r8), dimension(nw) :: fw =                                 &
     &            (/                                                    &
     &             0.25, 0.25, 0.25, 0.25                               &
     &              /)
# endif

      CONTAINS

      SUBROUTINE allocate_filter (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
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
      IF (ng.eq.1 ) allocate ( FILTER(Ngrids) )
!
      allocate ( FILTER(ng) % filzeta(LBi:UBi,LBj:UBj,FIL) )
      allocate ( FILTER(ng) % filu2d(LBi:UBi,LBj:UBj,FIL) )
      allocate ( FILTER(ng) % filv2d(LBi:UBi,LBj:UBj,FIL) )
# ifdef SOLVE3D
      allocate ( FILTER(ng) % filu3d(LBi:UBi,LBj:UBj,N(ng),FIL) )
      allocate ( FILTER(ng) % filv3d(LBi:UBi,LBj:UBj,N(ng),FIL) )
#  ifndef FILTRIM
      allocate ( FILTER(ng) % filrho(LBi:UBi,LBj:UBj,N(ng),FIL) )
      allocate ( FILTER(ng) % filt(LBi:UBi,LBj:UBj,N(ng),NAT,FIL) )
      allocate ( FILTER(ng) % filw3d(LBi:UBi,LBj:UBj,0:N(ng),FIL) )
#  endif
#  ifdef LMD_SKPP
      allocate ( FILTER(ng) % filhsbl(LBi:UBi,LBj:UBj,FIL) )
#  endif
#  ifdef LMD_BKPP
      allocate ( FILTER(ng) % filhbbl(LBi:UBi,LBj:UBj,FIL) )
#  endif
#  ifdef AVERAGES_AKV
      allocate ( FILTER(ng) % filAKv(LBi:UBi,LBj:UBj,0:N(ng),FIL) )
#  endif
#  ifdef AVERAGES_AKT
      allocate ( FILTER(ng) % filAKt(LBi:UBi,LBj:UBj,0:N(ng),FIL) )
#  endif
#  ifdef AVERAGES_AKS
      allocate ( FILTER(ng) % filAKs(LBi:UBi,LBj:UBj,0:N(ng),FIL) )
#  endif
#  ifdef AVERAGES_FLUXES
      allocate ( FILTER(ng) % filstf(LBi:UBi,LBj:UBj,FIL) )
      allocate ( FILTER(ng) % filswf(LBi:UBi,LBj:UBj,FIL) )
#   ifdef BULK_FLUXES
      allocate ( FILTER(ng) % fillhf(LBi:UBi,LBj:UBj,FIL) )
      allocate ( FILTER(ng) % filshf(LBi:UBi,LBj:UBj,FIL) )
      allocate ( FILTER(ng) % fillrf(LBi:UBi,LBj:UBj,FIL) )
#   endif
#   ifdef SHORTWAVE
      allocate ( FILTER(ng) % filsrf(LBi:UBi,LBj:UBj,FIL) )
#   endif
#  endif
#  ifdef ICE_MODEL
          allocate ( FILTER(NG) % filui(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filvi(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filai(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filhi(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filhsn(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filsfwat(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filtis(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filti(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filAgeice(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % fils0mk(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filt0mk(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filWfr(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filWai(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filWao(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filWio(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filWro(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filsig11(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filsig22(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filsig12(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filio_mflux(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filutau_iw(LBi:UBi,LBj:UBj,FIL) )
          allocate ( FILTER(NG) % filchu_iw(LBi:UBi,LBj:UBj,FIL) )
#  endif
# endif
# ifdef AVERAGES_FLUXES
      allocate ( FILTER(ng) % filsus(LBi:UBi,LBj:UBj,FIL) )
      allocate ( FILTER(ng) % filsvs(LBi:UBi,LBj:UBj,FIL) )
# endif

      RETURN
      END SUBROUTINE allocate_filter

      SUBROUTINE initialize_filter (ng, tile)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
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
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: i, j, if, k, itrc
      integer :: Imin, Imax, Jmin, Jmax
      real(r8), parameter :: IniVal = 0.0_r8
!
# include "set_bounds.h"

#ifdef _OPENMP 
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
#else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO if=1,FIL
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            FILTER(ng) % filzeta(i,j,if) = IniVal
            FILTER(ng) % filu2d(i,j,if) = IniVal
            FILTER(ng) % filv2d(i,j,if) = IniVal
# ifdef SOLVE3D
#  ifdef LMD_BKPP
            FILTER(ng) % filhbbl(i,j,if) = IniVal
#  endif
#  ifdef LMD_SKPP
            FILTER(ng) % filhsbl(i,j,if) = IniVal
#  endif
#  ifdef AVERAGES_FLUXES
            FILTER(ng) % filstf(i,j,if) = IniVal
            FILTER(ng) % filswf(i,j,if) = IniVal
#   ifdef BULK_FLUXES
            FILTER(ng) % fillhf(i,j,if) = IniVal
            FILTER(ng) % filshf(i,j,if) = IniVal
            FILTER(ng) % fillrf(i,j,if) = IniVal
#   endif
#   ifdef SHORTWAVE
            FILTER(ng) % filsrf(i,j,if) = IniVal
#   endif
#  endif
#  ifdef ICE_MODEL
            FILTER(ng) % filui(i,j,if) = IniVal
            FILTER(ng) % filvi(i,j,if) = IniVal
            FILTER(ng) % filai(i,j,if) = IniVal
            FILTER(ng) % filhi(i,j,if) = IniVal
            FILTER(ng) % filhsn(i,j,if) = IniVal
            FILTER(ng) % filsfwat(i,j,if) = IniVal
            FILTER(ng) % filtis(i,j,if) = IniVal
            FILTER(ng) % filti(i,j,if) = IniVal
            FILTER(ng) % filAgeice(i,j,if) = IniVal
            FILTER(ng) % fils0mk(i,j,if) = IniVal
            FILTER(ng) % filt0mk(i,j,if) = IniVal
            FILTER(ng) % filWfr(i,j,if) = IniVal
            FILTER(ng) % filWai(i,j,if) = IniVal
            FILTER(ng) % filWao(i,j,if) = IniVal
            FILTER(ng) % filWio(i,j,if) = IniVal
            FILTER(ng) % filWro(i,j,if) = IniVal
            FILTER(ng) % filsig11(i,j,if) = IniVal
            FILTER(ng) % filsig22(i,j,if) = IniVal
            FILTER(ng) % filsig12(i,j,if) = IniVal
            FILTER(ng) % filio_mflux(i,j,if) = IniVal
            FILTER(ng) % filutau_iw(i,j,if) = IniVal
            FILTER(ng) % filchu_iw(i,j,if) = IniVal
#  endif
# endif
# ifdef AVERAGES_FLUXES
            FILTER(ng) % filsus(i,j,if) = IniVal
            FILTER(ng) % filsvs(i,j,if) = IniVal
# endif
          END DO
        END DO
      END DO
# ifdef SOLVE3D
      DO if=1,FIL
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filu3d(i,j,k,if) = IniVal
              FILTER(ng) % filv3d(i,j,k,if) = IniVal
#  ifndef FILTRIM
              FILTER(ng) % filrho(i,j,k,if) = IniVal
#  endif
            END DO
            DO itrc=1,NAT
              DO i=Imin,Imax
                FILTER(ng) % filt(i,j,k,itrc,if)= IniVal
              END DO
            END DO
          END DO
        END DO
      END DO
      DO if=1,FIL
        DO k=0,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              FILTER(ng) % filw3d(i,j,k,if) = IniVal
#  ifdef AVERAGES_AKV
              FILTER(ng) % filAKv(i,j,k,if) = IniVal
#  endif
#  ifdef AVERAGES_AKT
              FILTER(ng) % filAKt(i,j,k,if) = IniVal
#  endif
#  ifdef AVERAGES_AKS
              FILTER(ng) % filAKs(i,j,k,if) = IniVal
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
      RETURN
      END SUBROUTINE initialize_filter
#endif
      END MODULE mod_filter

